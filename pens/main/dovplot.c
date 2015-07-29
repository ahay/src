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
 *  source file:   ./filters/dovplot.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */
/*
 * Joe Dellinger Oct 18 1987
 * 	Made text fatness scale with txscale and scale for GKS compatibility
 * Joe Dellinger Nov 9 1987
 *	Changed "\n" to CRLF in MESG_TEXT calls, since newlines are
 *	no longer mapped into Carriage-return linefeeds automatically.
 *	(This was due to Stew's turning off of output translation.)
 * Joe Dellinger Dec 9 1987
 *	Made the "bad color level" error message more verbose.
 * Joe Dellinger Dec 19 1987
 *	Call dev.attributes(NEW_DASH,...) when the dash pattern has
 *	been changed, call dev.attributes(NEW_PAT,...) when a new raster
 *	pattern has been loaded, call dev.attributes(NEW_FONT,...)
 *	when the txfont, txprec, or txovly has been changed,
 *	dev.attributes(NEW_OVERLAY,...) when the overlay mode has
 *	been changed, dev.attributes(NEW_ALIGN) when the text
 *	alignment mode has been changed, dev.attributes(NEW_FAT) when
 *	the fatness has been changed. Call dev.message(MESG_MESSAGE)
 *	when a message is generated via the VP_MESSAGE vplot command.
 * Joe Dellinger Dec 20 1987
 *	It's up to the device to turn around and call gentext if
 *	txfont < NUMGENFONT.
 * Joe Dellinger Jan 8 1988
 *	Added reset_parameters() to frontend and dovplot.
 *	Turn off ALL polygon shading if shade=NO, even that generated
 *	by raster, text, etc.
 * Joe Dellinger Jan 10 1988
 *	Don't blindly scale up fatmult each time dovplot thinks it's
 *	being called for the first time.
 * Joe Dellinger Jan 14 1988
 *	The text alignment should also be reset between plots.
 * Joe Dellinger Jan 20 1988
 *	Added VP_BEGIN_GROUP and VP_END_GROUP. Fixed nplots bug.
 * Joe Dellinger Jan 22 1988
 *	Added xdim_orig and ydim_orig to pat.h.
 * Joe Dellinger Jan 26 1988
 *	Fixed EOF bug in getvpstring(). Changed rule for mapping
 *	colors to grey levels to match how B+W TV's do it.
 * Joe Dellinger Feb 10 1988
 *	Weight color mismatches to tally with grey mapping method.
 * Joe Dellinger Feb 12 1988
 *	Make sure that VP_WINDOW commands that specify an inverted
 *	clipping rectangle cause everything to be clipped away.
 * Joe Dellinger Feb 16 1988
 *	Make number of arguments passed to dev.attributes and dev.raster
 *	constant to head off future catastrophe from Sun IV's.
 * Joe Dellinger Feb 22 1988
 *	Created INT_PAUSE to be separate from INT_GET_STRING.
 * Joe Dellinger Feb 24 1988
 *	Txfont, txprec, txovly, style, dash line style also reset
 *	on erases.
 *	Only the erase command is handled specially at the start
 *	of the file. This was causing initial setstyle commands
 *	to be left out of their groups!
 * Joe Dellinger Feb 26 1988
 *	Allow vpattributes(END_GROUP,...) to signal dovplot whether it
 *	should exit or not after the call.
 * Joe Dellinger Feb 28 1988
 *	Fatness behaves like a geometric attribute.
 *	Can't use the *= operator on the Sun!
 *	if integer i = 10, then i *= .5 sets i to ZERO.
 *	I wish Sun would get their act together.
 *	This explains several mysterious bugs on the Suns.
 * Joe Dellinger Mar 17 1988
 *	Added "colormask" option.
 * Joe Dellinger March 28 1988
 *	Code to be efficient and turn off dithering when the input data
 *	happens to already be black-and-white introduces a bug.
 *	Best solution is to get rid of the special case. While less
 *	efficient computer-wise, it is much better code-wise.
 * Joe Dellinger March 29 1988
 *	Draw the polygon border AFTER you draw the polygon, not before.
 * Joe Dellinger April 15 1988
 *	Make colormask work more reasonably.
 * Joe Dellinger April 17 1988
 *	To be consistent between Color and Monochrome devices,
 *	DEFAULT_COLOR must be COLOR_MAPped just like any other
 *	color, internally generated or externally. Fixed bug
 *	in COLOR_MAP macro definition (; on the end! oops!)
 * Joe Dellinger April 29 1988
 *	Txscale should be applied consistently for all 3 kinds
 *	of text.
 * Joe Dellinger Oct 28 1988
 *	Put in "txsquare=no" option.
 * Joe Dellinger May 7 1989
 * 	Put in break=ignore option.
 * Joe Dellinger August 2 1989
 *	Put in redmap, greenmap, bluemap, redpow, greenpow, bluepow
 *	options.
 * Dave Nichols  Feb. 21 1991
 *	Add extern declaration of geth (it returns a short).
 * Dave Nichols  Apr. 12 1991
 *	Add check on return value of dev.interact. See if it is DOVPLOT_EXIT.
 * Joe Dellinger Apr. 2 1992
 *	Added argument to dev.attributes(END_GROUP,...) to give more info
 *	to the device so dovplot can be used a frame renderer using the
 * 	begin and endgroup calls to keep track of where each frame begins
 *	and ends.
 * Dave Nichols 9-17-92
 *      Include stdlib.h if it is available instead of defining malloc.
 * Joe Dellinger, David Lumley 10-03-94
 *	Created the external variable "ras_allgrey" to allow smart color
 *	devices (like postscript) to know whether a color raster will
 *      actually require color to plot.
 */

/*
 * process a vplot file or stream
 * Keywords: vplot pen generic
 */

#include 	<stdio.h>
#include 	<math.h>
#include	<ctype.h>
#include	<string.h>

#include <unistd.h>

#include <rsfplot.h>

#include	"../include/params.h"	/* for machine dependencies */
#include	"../include/enum.h"
#include	"../include/err.h"
#include	"../include/attrcom.h"
#include	"../include/intcom.h"
#include	"../include/mesgcom.h"
#include	"../include/erasecom.h"
#include	"../include/closestat.h"
#include	"../include/getxy.h"
#include	"../include/pat.h"
#include	"../include/vertex.h"
#include	"../include/extern.h"
#include	"../include/round.h"
#include	"../include/device.h"

#include "../genlib/genpen.h" 
#include "../utilities/util.h"

#include "dovplot.h"
#include "init_vplot.h"

#define COLOR_MAP(A) (color_set[(A)][MAP])
#define GREY_MAP(A) (color_set[(A)][_GREY])

/*
 * Banished to an include file so that "indent" can't get its hands on it
 */
#include "../include/readraster.h"

void end_of_file(void);
void set_table(void);
extern bool      wantras, smart_raster;
extern int      brake;
extern FILE    *controltty;
extern int      cur_color;
extern int      pat_color;
extern int      epause;
extern int      erase;
extern int      ever_called;
extern int      fatbase;
extern int      ifat;
extern int      first_time;
extern bool      framewindows;
extern int      next_color;
extern int      ipat;
extern int      nplots;
extern bool      overlay;
extern struct pat pat[];
extern char     group_name[];
extern int      group_number;
extern char    *txbuffer;
extern int      txbuflen;
extern struct vertex *vxbuffer;
extern int      vxbuflen;
extern int      window;
extern int      xwmax, xwmin, ywmax, ywmin;
extern int      xWmax, xWmin, yWmax, yWmin;
int             xwmin_last, xwmax_last, ywmin_last, ywmax_last;
extern int      xnew, ynew;
extern int      xold, yold;
extern int      xorigin, yorigin;
extern float    scale;
extern float    xscale;
extern float    yscale;
extern vp_plotstyle default_style;
extern int      default_txfont, default_txprec, default_txovly;
extern bool     default_overlay;
extern int      color_set[MAX_COL + 1][_NUM_PRIM];
extern int      greycorr ();
extern int      num_col_8;
extern int      xret, yret;
extern char     interact[];
extern float    dashsum;
extern float    dashpos;
extern float    dashes[];
extern bool     colormask[];
extern float    redmap[4], greenmap[4], bluemap[4];
extern float    redpow, greenpow, bluepow;
int             ras_allgrey = YES;

#include <stdlib.h>

long int        ftell ();
extern void     wlimit ();

int             need_devcolor = NO;

static void reset (void);
static void outline_window (void);
static void update_color (void);
static void getvpstring (void);
static int getpolygon (int npts);
static void fill_background (void);

void gen_dovplot (int nn, FILE **inpltin, char** innames)
/*< Call dovplot, cycling through the input plot files >*/
{
    int             ii;

    for (ii = 0; ii < nn; ii++)
    {
/*
 * Set things up for dovplot.
 * The "-5" is to leave room for dovplot to add a little onto the
 * end of pltname if it wants.
 */
	pltin = inpltin[ii];
	strncpy (pltname, innames[ii], MAXFLEN - 5);
	dovplot ();
	fclose (pltin);
    }
}

void dovplot (void)
/*< do vplot >*/
{
    int             i, j, k, c;
    int             key, size, npts;
    int             nmul;
    int             nx, ny, nx_orig, ny_orig;
    int             orient, ras_orient;
    register int   *ptr=NULL;
    int             nx_mult, ny_mult;
    int             nx_temp, ny_temp;
    int            *tempbuf, *ptemp;
    vp_plotstyle    new_style;
    int             starterase = 0;
    int             hacol[NHATCH * 2], hafat[NHATCH * 2], haoff[NHATCH * 2],
	hasiz[NHATCH * 2], numhatch;
    float           angle, xrasmult, yrasmult;
    int             xpix, ypix, num_pat, num_byte;
    int             yrast, lastrast;
    int             xvr_min, xvr_max, yvr_min, yvr_max;
    int             xr_min, xr_max, yr_min, yr_max;
    int             xvru_min, xvru_max, yvru_min, yvru_max;
    int             xru_min, yru_max;
    int             xxx[4], yyy[4];
    int             pos, ii, jj, kk, num_rep, ras_offset, dither_it;
    unsigned char  *rasterline, *rasterline2, **outraster, **outraster2;
    unsigned char   ibyte;
    int             xnewer, ynewer;
    int             xtext0, xtext1, xtext2, ytext0, ytext1, ytext2;
    int             type;
    int            *marker_vec, *mvec;
    int             savefat;
    float           savefatmult;
    char            string[MAXFLEN + 1];

    /*
     * Check to make sure we really got anything before we claim we've made a
     * plot.
     */
    if ((c = getc (pltin)) == EOF)
	return;
    ungetc ((char) c, pltin);


    while (((c = getc (pltin)) == VP_ERASE) || 
	   ((c == VP_BREAK) && (dev.brake == BREAK_ERASE)))
    {
	/*
	 * This is a baby version of the main switch that takes up most of
	 * this file.
	 */
	switch (c)
	{
	    case VP_ERASE:
	    case VP_BREAK:
		starterase = 1;
		break;
	}
    }

    /*
     * Files that consist only of erase characters and nothing else are
     * ignored.
     */
    if (c == EOF)
	return;
    ungetc ((char) c, pltin);

    if (first_time)
    {
	dev.reset ();
/*
 * Device is now officially open for orders.
 */
	if (dev.xmax <= dev.xmin ||
	    dev.ymax <= dev.ymin ||
	    dev.pixels_per_inch == 0. ||
	    dev.aspect_ratio == 0. ||
	    dev.num_col == -1)
	    ERR (FATAL, name, "Critical variables left unset by device!");

/*
 * Set maximum clipping window for dev.attributes(SET_WINDOW)
 */
	xwmax_last = dev.xmax;
	xwmin_last = dev.xmin;
	ywmax_last = dev.ymax;
	ywmin_last = dev.ymin;

/*
 * Set up color maps
 */
	init_colors ();

	ever_called = 1;
    }
    else
    {
	dev.close (CLOSE_FLUSH);
	if (epause > 0)
	{
	    sleep ((unsigned) epause);
	}
	else if (epause < 0)
	{
	    dev.close (CLOSE_FLUSH);
	    message (MESG_ERASE,NULL);
	    message (MESG_ON,NULL);
	    message (MESG_HOME,NULL);
	    message (MESG_READY,NULL);
	    message (MESG_HIGHLIGHT_ON,NULL);
	    message (MESG_TEXT, "Type Return to Continue...  ");
	    message (MESG_DONE,NULL);
/*	    dev.interact (INT_PAUSE, controltty, string);*/
	    ii = dev.interact (INT_PAUSE, controltty, string);
	    if (ii == DOVPLOT_EXIT)
	    {
		dev.close (CLOSE_FLUSH);
		return;
	    }
	    message (MESG_HIGHLIGHT_OFF,NULL);
	    message (MESG_DONE,NULL);
	    message (MESG_OFF,NULL);
	    message (MESG_ERASE,NULL);
	}
	/*
	 * Inquire point back from device
	 */
	if (interact[0] != '\0')
	{
	    getapoint ();
	}
    }

    /*
     * Erase the screen to background color
     */
    if (((erase & FORCE_INITIAL) && (first_time || (erase & DO_LITERALS)))
	|| (starterase && (erase & DO_LITERALS)))
    {

	if (first_time)
	    dev.erase (ERASE_START);
	else
	{
	    dev.erase (ERASE_MIDDLE);
	    nplots++;
	}
	dev.close (CLOSE_FLUSH);
    }

    first_time = NO;

/*
 * Reset fatness, cur_color, etc.
 */
    new_style = default_style;
    setstyle (new_style);
    reset ();

/*
 * Make SURE the color is what it's supposed to be, just to be safe.
 */
    dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
    need_devcolor = NO;

    message (MESG_OFF,NULL);
    dev.close (CLOSE_FLUSH);

/*
 * Start a group that will contain this frame (always group 0).
 */
    ii = ftell (pltin);
    sprintf (group_name, "%s.%d", pltname, nplots);
    dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
    group_number++;

    if (framewindows)
	outline_window ();

/*
 * Finally, here's the main loop that does the actual processing of vplot.
 * Notice that some tricky commands (VP_DRAW, VP_SET_COLOR_TABLE) look ahead
 * at the next command to try to anticipate often-occuring situations,
 * so things get a little complicated here.
 */
    while ((c = getc (pltin)) != EOF)
    {
	switch (c)		/* command list */
	{
	    case VP_SETSTYLE:	/* set the style */
		c = getc (pltin);
		if ((c == 'r') || (c == 'R') || (c == 'm') || (c == 'M'))
		    new_style = VP_ROTATED;
		else if ((c == 'o') || (c == 'O'))
		    new_style = VP_OLD;
		else if ((c == 'a') || (c == 'A'))
		    new_style = VP_ABSOLUTE;
		else
		    new_style = VP_STANDARD;
		setstyle (new_style);

		if (framewindows)
		    outline_window ();

		break;
	    case VP_MOVE:		/* move */
		/*
		 * Reset position in dash pattern.
		 */
		dashpos = 0.;

		GETXY (xold, yold);
		break;
	    case VP_DRAW:		/* draw */
		GETXY (xnew, ynew);
		update_color ();
		while ((c = getc (pltin)) == VP_DRAW)
		{
		    GETXY (xnewer, ynewer);
		    /*
		     * Is it the same point?
		     */
		    if (xnewer == xnew && ynewer == ynew)
			continue;
		    /*
		     * Is it colinear and horizontal or vertical?
		     */
		    if ((ynewer == ynew && ynew == yold &&
			 ((xnewer > xnew) == (xnew > xold))) ||
			(xnewer == xnew && xnew == xold &&
			 ((ynewer > ynew) == (ynew > yold))))
		    {
			ynew = ynewer;
			xnew = xnewer;
			continue;
		    }
		    dev.vector (xold, yold, xnew, ynew, fat, dashon);
		    xold = xnew;
		    yold = ynew;
		    xnew = xnewer;
		    ynew = ynewer;
		}
		dev.vector (xold, yold, xnew, ynew, fat, dashon);
		xold = xnew;
		yold = ynew;
		if (c == EOF)
		    goto End_of_file;
		ungetc ((char) c, pltin);
		break;
	    case VP_PLINE:		/* polyline */
		/*
		 * Reset position in dash pattern.
		 */
		dashpos = 0.;

		npts = geth (pltin);
		if (npts == 0)
		    break;
		GETXY (xold, yold);
		npts--;
		if (npts == 0)
		    break;
		GETXY (xnew, ynew);
		npts--;
		update_color ();
		while (npts > 0)
		{
		    GETXY (xnewer, ynewer);
		    npts--;
		    /*
		     * Is it the same point?
		     */
		    if (xnewer == xnew && ynewer == ynew)
			continue;
		    /*
		     * Is it colinear and horizontal or vertical?
		     */
		    if ((ynewer == ynew && ynew == yold &&
			 ((xnewer > xnew) == (xnew > xold))) ||
			(xnewer == xnew && xnew == xold &&
			 ((ynewer > ynew) == (ynew > yold))))
		    {
			ynew = ynewer;
			xnew = xnewer;
			continue;
		    }
		    dev.vector (xold, yold, xnew, ynew, fat, dashon);
		    xold = xnew;
		    yold = ynew;
		    xnew = xnewer;
		    ynew = ynewer;
		}
		dev.vector (xold, yold, xnew, ynew, fat, dashon);
		xold = xnew;
		yold = ynew;
		break;
	    case VP_PMARK:		/* polymarker */
		npts = geth (pltin);/* how many markers ? */
		type = geth (pltin);/* what symbol? (any positive integer) */
		size = geth (pltin);/* How big? */
		size = size * mkscale * vdevscale * RPERIN / TXPERIN;

		if (npts == 0)
		    break;
		/* allocate space for the points */
		marker_vec = sf_intalloc (npts * 2);
		mvec = marker_vec;
		if (mvec == NULL)
		    ERR (FATAL, name, "Can't malloc memory for markers!");

		/* read the locations, and transform them to device coordinates */
		for (ii = 0; ii < npts; ii++)
		{
		    GETXY (xnewer, ynewer);
		    *mvec = xnewer;
		    ++mvec;
		    *mvec = ynewer;
		    ++mvec;
		}
		update_color ();
		/* call the device routine to display the markers */
		dev.marker (npts, type, size, marker_vec);
		/* release the storage used for the points */
		free ((char *) marker_vec);
		break;
	    case VP_ORIGIN:	/* set origin */
		dev.xorigin = geth (pltin);
		dev.yorigin = geth (pltin);
		break;
	    case VP_BEGIN_GROUP:
		ii = ftell (pltin) - 1;
		getvpstring ();
		strncpy (group_name, txbuffer, MAXFLEN);
		group_name[MAXFLEN] = '\0';
		dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
		group_number++;
		break;
	    case VP_END_GROUP:
		group_number--;
		if (group_number < 1)
		{
		    ERR (WARN, name,
			 "group invalidly nested");
		    group_number = 1;
		}
		else
		{
		    dev.attributes (END_GROUP, group_number, USER_EGROUP, 0, 0);
		}
		break;
	    case VP_GTEXT:		/* GKS-like text */
		xtext0 = 0;
		ytext0 = 0;
		vptodevxy_text (xtext0, ytext0, &xtext0, &ytext0);
		GETXY_TEXT (xtext1, ytext1);
		GETXY_TEXT (xtext2, ytext2);
	    g_text:
		xtext1 -= xtext0;
		xtext2 -= xtext0;
		ytext1 -= ytext0;
		ytext2 -= ytext0;

		savefat = fat;
		savefatmult = fatmult;
		fatmult *= txscale;
		if (ifat >= 0)
		{
		    fat = fatmult * (float) (ifat + fatbase);
		}
		else
		{
		    fat = -1;
		}
		update_color ();
		getvpstring ();
/*
 * Fonts less than NUMGENFONT reserved for gentext fonts:
 * up to the device to enforce that rule, though.
 */
		dev.text (txbuffer,
			  txscale * (float) xtext1 / TEXTVECSCALE,
			  txscale * (float) ytext1 / TEXTVECSCALE,
			  txscale * (float) xtext2 / TEXTVECSCALE,
			  txscale * (float) ytext2 / TEXTVECSCALE);
		fat = savefat;
		fatmult = savefatmult;
		break;
	    case VP_TEXT:		/* text */
		size = geth (pltin);
		size = size * (float) RPERIN / (float) TXPERIN;
		orient = (int) geth (pltin);
	    new_text:
		xtext0 = 0;
		ytext0 = 0;
		/* Character path direction */
		xtext1 =
		    ROUND (TEXTVECSCALE * size * cos (orient * 3.14159 / 180.));
		ytext1 =
		    ROUND (TEXTVECSCALE * size * sin (orient * 3.14159 / 180.));
		/* Character up vector direction */
		orient += 90;
		xtext2 =
		    ROUND (TEXTVECSCALE * size * cos (orient * 3.14159 / 180.));
		ytext2 =
		    ROUND (TEXTVECSCALE * size * sin (orient * 3.14159 / 180.));
		vptodevxy_text (xtext0, ytext0, &xtext0, &ytext0);
		vptodevxy_text (xtext1, ytext1, &xtext1, &ytext1);
		vptodevxy_text (xtext2, ytext2, &xtext2, &ytext2);
		goto g_text;
		break;
	    case VP_OLDTEXT:	/* archaic format text */
		if ((key = geth (pltin)) < 0)
		    ERR (FATAL, name, "invalid text key");
		size = (key & 037);
		size = size * (float) RPERIN / (float) TXPERIN;
		orient = (int) (((key & 0140) >> 5) * 90);
		goto new_text;
		break;
	    case VP_OLDAREA:	/* polygon */
		npts = geth (pltin);
		afat = geth (pltin);
		if (afat >= 0)
		{
		    afat = fatmult * (fatbase + afat);
		}
		nx_temp = geth (pltin);
		ny_temp = geth (pltin);

/*
 * If a monochrome device, then fill with dot pattern.
 * If a color device, fill with solid color.
 */
		nx = nx_temp * patternmult;
		if (nx_temp == 1 || (nx_temp > 0 && nx == 0))
		    nx = 1;
		ny = ny_temp * patternmult;
		if (ny_temp == 1 || (ny_temp > 0 && ny == 0))
		    ny = 1;

		nx_orig = nx;
		ny_orig = ny;

		if (!mono)
		{
		    if (nx_temp == 0 || ny_temp == 0)
		    {
			nx = 0;
			ny = 0;
		    }
		    else
		    {
			nx = 1;
			ny = 1;
		    }
		}
		/*
		 * Create a temporary pattern
		 */
		ipat = 0;
		if (nx * ny > 0)
		{
		    if ((ptr = (int *) calloc ((unsigned) nx * ny, sizeof (int))) == NULL)
			ERR (FATAL, name, "cannot alloc memory to load pattern");
		    pat[ipat].patbits = ptr;
		    ptr[(nx * ny) - 1] = cur_color;
		}
		else
		    pat[ipat].patbits = NULL;

		pat[ipat].xdim = nx;
		pat[ipat].ydim = ny;
		pat[ipat].xdim_orig = nx_orig;
		pat[ipat].ydim_orig = ny_orig;
		npts = getpolygon (npts);
		update_color ();

/* Do the polygon */
		if (npts > 2 && shade && nx > 0 && ny > 0)
		    dev.area (npts, vxbuffer);

/* And then its border */
		if (afat >= 0)
		    vecoutline (vxbuffer);

		if (nx * ny > 0)
		    free ((char *) ptr);
		break;
	    case VP_AREA:		/* polygon fill */
		npts = geth (pltin);
		ipat = pat_color;

		if (pat[ipat].patbits == NULL && shade)
		{
		    /*
		     * Create a default pattern (solid fill with this color) If
		     * you don't like this default for a particular device, ie, a
		     * black and white one, then the device itself can set up
		     * defaults for colors 0 through 7 in dev.open
		     */
		    nx = 1;
		    ny = 1;

		    if ((ptr = (int *) calloc ((unsigned) nx * ny, sizeof (int))) == NULL)
			ERR (FATAL, name, "cannot alloc memory to load pattern");
		    pat[ipat].patbits = ptr;
		    ptr[(nx * ny) - 1] = cur_color;
		    pat[ipat].xdim = nx;
		    pat[ipat].ydim = ny;
		    pat[ipat].xdim_orig = nx;
		    pat[ipat].ydim_orig = ny;
		}

		npts = getpolygon (npts);

		if (!shade)
		{
		    /* At least draw the boundary to show where it is */
		    afat = 0;
		    update_color ();
		    vecoutline (vxbuffer);
		}
		else
		{
		    /* See whether raster or hatch area */
		    if (pat[ipat].xdim >= 0)
		    {
			/* raster */
			if (npts > 2 && pat[ipat].ydim > 0 && pat[ipat].xdim > 0)
			{
			    update_color ();
			    dev.area (npts, vxbuffer);
			}
		    }
		    else
		    {
			/* hatch */
			numhatch = -pat[ipat].xdim;
			angle = 3.14159 * pat[ipat].ydim / 180.;
			ptr = pat[ipat].patbits;
			for (i = 0; i < numhatch * 2; i++)
			{
			    hafat[i] = (*ptr++);
			    hacol[i] = (*ptr++);
			    haoff[i] = (*ptr++);
			    hasiz[i] = (*ptr++);
			}
			if (npts > 2)
			{
			    /*
			     * Hatch patterns don't rotate. The polygon does, the
			     * fill pattern doesn't.
			     */
			    genhatch (npts, numhatch, angle, hafat, hacol, haoff, hasiz, vxbuffer);
			}
		    }
		}
		break;
	    case VP_FAT:		/* fat */
		/*
		 * negative fat always means to not draw the line at all.
		 */
		ifat = geth (pltin);
		if (ifat >= 0)
		{
		    fat = fatmult * (float) (ifat + fatbase);
		}
		else
		{
		    fat = -1;
		}
		dev.attributes (NEW_FAT, fat, 0, 0, 0);
		break;
	    case VP_COLOR:		/* change color */
		pat_color = geth (pltin);
		if (pat_color > MAX_COL || pat_color < 0)
		{
		    ERR (WARN, name, "bad color number %d (max %d, min 0)",
			 pat_color, MAX_COL);
		    pat_color = DEFAULT_COLOR;
		}
		next_color = COLOR_MAP (pat_color);
		if (next_color != cur_color)
		{
		    cur_color = next_color;
		    need_devcolor = YES;
		}
		/*
		 * Pattern zero reserved for the OLD_AREA command, so increment
		 * the rest by one to make room.
		 */
		pat_color++;
		break;
	    case VP_SET_COLOR_TABLE:	/* set color table entry */
		set_table();
		break;
	    case VP_PURGE:		/* purge pltout buffers */
		dev.close (CLOSE_FLUSH);
		break;
	    case VP_BACKGROUND:
		if (honor_background) {
		    if (dev.smart_background) {
			dev.erase (ERASE_BACKGROUND);
		    } else {
			fill_background();
		    }
		}
		break;
	    case VP_BREAK:		/* break */
		if (dev.brake == BREAK_IGNORE)
		    break;
		/* NOTE break IS IN AN "if": WE USUALLY FALL THROUGH HERE! */
	    case VP_ERASE:		/* erase (and break falls through here too) */
		dev.close (CLOSE_FLUSH);

		/*
		 * Erases and breaks can't occur inside groups.
		 */
		group_number--;
		if (group_number != 0)
		{
		    ERR (WARN, name,
			 "group contains erase or break");
		    group_number = 0;
		}
		else
		{
		    if (erase & DO_LITERALS)
		    {
			if (c == VP_ERASE || dev.brake)
			    dev.attributes (END_GROUP, group_number, ERASE_EGROUP, 0, 0);
			else
			    dev.attributes (END_GROUP, group_number, BREAK_EGROUP, 0, 0);
		    }
		    else
			dev.attributes (END_GROUP, group_number, IGNORE_EGROUP, 0, 0);
		}

		if (epause < 0)
		{
		    message (MESG_ERASE,NULL);
		    message (MESG_ON,NULL);
		    message (MESG_HOME,NULL);
		    message (MESG_READY,NULL);
		    message (MESG_HIGHLIGHT_ON,NULL);
		    message (MESG_TEXT, "Type Return to Continue...  ");
		    message (MESG_DONE,NULL);
		    ii = dev.interact (INT_PAUSE, controltty, string);
		    if (ii == DOVPLOT_EXIT)
			return;
		    message (MESG_HIGHLIGHT_OFF,NULL);
		    message (MESG_DONE,NULL);
		    message (MESG_OFF,NULL);
		    message (MESG_ERASE,NULL);
		}
		else
		{
		    if (epause > 0)
			sleep ((unsigned) epause);
		}
		/*
		 * Inquire point back from device
		 */
		if (interact[0] != '\0')
		{
		    getapoint ();
		}

		if (erase & DO_LITERALS) {
		    if ((c == VP_ERASE) || dev.brake)
		    {
			dev.erase (ERASE_MIDDLE);
			nplots++;
		    }
		    else
		    {
			dev.erase (ERASE_BREAK);
			nplots++;
		    }
		}

		new_style = default_style;
		setstyle (new_style);
		reset ();

		/*
		 * Start a new group level 0 to contain the next frame (separated
		 * from the previous one by either an erase or break)
		 */
		ii = ftell (pltin);
		sprintf (group_name, "%s.%d", pltname, nplots);
		dev.attributes (BEGIN_GROUP, group_number, ii, 0, 0);
		group_number++;

		if (framewindows)
		    outline_window ();

		break;
	    case VP_WINDOW:	/* window */
		if (window)
		{
		    xwmin = geth (pltin);
		    ywmin = geth (pltin);
		    xwmax = geth (pltin);
		    ywmax = geth (pltin);

		    if (xwmin > xwmax || ywmin > ywmax)
		    {
/*
 * vptodevw will always return a window that is the "right way around",
 * even if the original window was inverted. So produce an inside-out
 * window which will clip everything away to nothing.
 */

			xwmin = xWmax + 1;
			xwmax = xWmin - 1;
			ywmin = yWmax + 1;
			ywmax = yWmin - 1;
		    }
		    else
		    {

			vptodevw (xwmin, ywmin, xwmax, ywmax, &xwmin, &ywmin, &xwmax, &ywmax);

			wlimit (xWmin, xWmax, &xwmin, &xwmax);
			wlimit (yWmin, yWmax, &ywmin, &ywmax);
		    }
		}
		else
		{
		    geth (pltin);
		    geth (pltin);
		    geth (pltin);
		    geth (pltin);
		}
		reset_windows ();
		if (framewindows)
		    outline_window ();
		break;
	    case VP_NOOP:		/* no op */
		break;
	    case VP_TXALIGN:	/* set text alignment */
		/*
		 * These are made available to the device text routine
		 */
		txalign.hor = geth (pltin);
		txalign.ver = geth (pltin);
		dev.attributes (NEW_ALIGN, txalign.hor, txalign.ver, 0, 0);
		break;
	    case VP_TXFONTPREC:	/* set text font */
		/*
		 * These are made available to the device text routine
		 */
		ii = geth (pltin);

		/*
		 * Force font substitution if serif fonts aren't OK
		 */
		if (! serifs_OK)
		{
		    if (ii == DEFAULT_HARDCOPY_FONT)
		    {
			ii = DEFAULT_SANSSERIF_FONT;
		    }
		    else if (ii == GREEK_SERIF_FONT)
		    {
			ii = GREEK_SANSSERIF_FONT;
		    }
		}

		if (ii >= 0)
		    dev.txfont = ii;
		else
		    ii = -1;

		jj = geth (pltin);
		if (jj >= 0)
		    dev.txprec = jj;
		else
		    jj = -1;

		kk = geth (pltin);
		if (kk >= 0)
		    dev.txovly = kk;
		else
		    kk = -1;

		/*
		 * Another way for the device to keep track of changes.
		 */
		dev.attributes (NEW_FONT, ii, jj, kk, 0);
		break;
	    case VP_OVERLAY:	/* change overlay mode */
		/*
		 * This is made available to the device dependent subroutines
		 */
		overlay = (bool) (1 == geth (pltin));
		dev.attributes (NEW_OVERLAY, overlay, 0, 0, 0);
		break;
	    case VP_PATLOAD:	/* load a pattern */
		nmul = geth (pltin);
		ny = geth (pltin);

		/* See whether Raster or Hatch */
		if (ny >= 0)
		{
		    /* Raster */
		    /* nmul gives pixels_per_inch pattern is designed for */
		    ny_mult = ny * (dev.pixels_per_inch / nmul) * patternmult / dev.aspect_ratio;
		    if (ny_mult == 0 && ny > 0)
			ny_mult = 1;
		    nx = geth (pltin);
		    nx_mult = nx * (dev.pixels_per_inch / nmul) * patternmult;
		    if (nx_mult == 0 && nx > 0)
			nx_mult = 1;

		    ipat = geth (pltin) + 1;
		    if (ipat > NPAT || ipat < 1)
			ERR (FATAL, name, "bad pattern number %d (max %d, min 1)", ipat, NPAT);
		    /*
		     * Free up pattern that may already be there
		     */
		    if (pat[ipat].patbits != NULL)
		    {
			free ((char *) pat[ipat].patbits);
		    }

		    if (nx_mult * ny_mult > 0)
		    {
			ptr = sf_intalloc (nx_mult * ny_mult);
			pat[ipat].patbits = ptr;
		    }
		    else
		    {
			ptr = sf_intalloc (1);
			pat[ipat].patbits = ptr;
			ptr[0] = 0;
		    }

		    if (nx * ny > 0)
		    {
			tempbuf = sf_intalloc (nx * ny);
		    }
		    else
			tempbuf = NULL;

		    /*
		     * read in pattern
		     */
		    ptemp = tempbuf;
		    for (j = 0; j < nx * ny; j++)
		    {
			k = geth (pltin);
			if (k > MAX_COL || k < 0)
			{
			    ERR (WARN, name, "bad color number in pattern %d (max %d, min 0)",
				 k, MAX_COL);
			    k = DEFAULT_COLOR;
			}
			*ptemp++ = COLOR_MAP (k);
		    }

		    /*
		     * copy into patbits, with "stretching"
		     */
		    for (i = 0; i < ny_mult; i++)
		    {
			ny_temp = i * ny / ny_mult;
			for (j = 0; j < nx_mult; j++)
			{
			    nx_temp = j * nx / nx_mult;
			    ptr[j + nx_mult * i] = tempbuf[nx_temp + nx * ny_temp];
			}
		    }
		    /*
		     * set dimensions of pattern
		     */
		    pat[ipat].xdim = nx_mult;
		    pat[ipat].ydim = ny_mult;
		    pat[ipat].xdim_orig = nx_mult;
		    pat[ipat].ydim_orig = ny_mult;

		    if (tempbuf != NULL)
			free ((char *) tempbuf);

		    dev.attributes (NEW_PAT, ipat, 0, 0, 0);
		}
		else
		{
		    /* Hatch Pattern */
		    /* nmul gives angle, ny is merely a flag */
		    nx = geth (pltin);
		    if (nx <= 0 || nx * 2 > NHATCH)
			ERR (FATAL, name, "bad numhatch %d (max %d/2, min 1)", nx, NHATCH);
		    ipat = geth (pltin) + 1;
		    if (ipat > NPAT || ipat < 1)
			ERR (FATAL, name, "bad pattern number %d (max %d, min 1)", ipat, NPAT);
		    /*
		     * Free up pattern that may already be there
		     */
		    if (pat[ipat].patbits != NULL)
		    {
			free ((char *) pat[ipat].patbits);
		    }
		    ptr = sf_intalloc ((unsigned) (nx * 4 * 2));
		    pat[ipat].patbits = ptr;

		    for (i = 0; i < nx * 2; i++)
		    {
			hafat[i] = geth (pltin);
			if (hafat[i] >= 0)
			{
			    hafat[i] = fatmult * (fatbase + hafat[i]);
			}
			k = geth (pltin);
			if (k > MAX_COL || k < 0)
			{
			    ERR (WARN, name, "bad color number in hatch %d (max %d, min 0)",
				 k, MAX_COL);
			    k = DEFAULT_COLOR;
			}
			hacol[i] = COLOR_MAP (k);
			haoff[i] = geth (pltin) * patternmult * dev.pixels_per_inch / RPERIN;
			hasiz[i] = geth (pltin);
		    }
		    /*
		     * Find the smallest hatch interval. 1/2 that, and then force
		     * everything to be a multiple of that. If this were not
		     * done, then hatch intervals that originally differed by
		     * some simple fraction might end up slightly off, causing
		     * very unsightly beats.
		     */
		    /*
		     * Upper quantization limit of 1/10 inch
		     */
		    k = (RPERIN / 10) * 2;
		    for (i = 0; i < nx * 2; i++)
		    {
			if (hasiz[i] > 0 && hasiz[i] < k)
			    k = hasiz[i];
		    }
		    j = k * (patternmult * dev.pixels_per_inch / RPERIN) / 2.;
		    if (j < 1)
			j = 1;
		    for (i = 0; i < nx * 2; i++)
		    {
			hasiz[i] = ((int) ((hasiz[i] * 2) / k)) * j;
		    }
		    /*
		     * The above algorithm also means that you can't have a hatch
		     * pattern come out on any device with a repetition rate of
		     * faster than once per two pixels.
		     */

		    for (i = 0; i < nx * 2; i++)
		    {
/*
 * Make sure haoff < hasiz
 */
			if (haoff[i] >= hasiz[i])
			    haoff[i] = 0;
			*ptr++ = hafat[i];
			*ptr++ = hacol[i];
			*ptr++ = haoff[i];
			*ptr++ = hasiz[i];
		    }
		    /*
		     * set numhatch and angle... dimensions are 2 * numhatch by 4
		     * . pat[ipat].xdim negative is a flag that this is a hatch
		     * pattern, not raster, and so must be treated differently.
		     */
		    /*
		     * numhatch, with neg as flag for hatch numhatch = 0 is OK,
		     * because that means don't fill, same as nx = 0.
		     */
		    pat[ipat].xdim = -nx;
		    pat[ipat].ydim = nmul;	/* angle */
		}
		break;
	    case VP_BIT_RASTER:	/* bit raster data */
	    case VP_BYTE_RASTER:	/* byte raster data */
		ras_orient = geth (pltin);
		if (rotate % 90 != 0)
		{
		    if (wantras)
		    {
			ERR (WARN, name, 
			     "Raster only possible in 4 principal orientations");
			wantras = false;
		    }
		}
		else
		{
		    ras_orient += rotate / 90;
		    if (ras_orient >= 0)
			ras_orient = ras_orient % 4;
		    else
			ras_orient = ((ras_orient % 4) + 4) % 4;
		}

		/*
		 * They're on their honor to not go out of bounds. This check is
		 * just for things that HAVE to go out of bounds.
		 */
		ras_offset = geth (pltin);

		if (ras_offset + 0 > MAX_COL || ras_offset + 255 < 0)
		{
		    ERR (FATAL, name, "Absurd raster offset %d", ras_offset);
		}

		xvr_min = geth (pltin);
		yvr_min = geth (pltin);
		xvr_max = geth (pltin);
		yvr_max = geth (pltin);
		vptodevw (xvr_min, yvr_min, xvr_max, yvr_max,
			  &xvr_min, &yvr_min, &xvr_max, &yvr_max);
		xvru_min = xvr_min;
		yvru_min = yvr_min;
		xvru_max = xvr_max;
		yvru_max = yvr_max;

		xpix = geth (pltin);
		ypix = geth (pltin);

		switch (ras_orient)
		{
		    case 0:
			xrasmult = (float) xpix / (float) (xvr_max - xvr_min);
			yrasmult = (float) ypix / (float) (yvr_max - yvr_min);
			xvr_max--;
			yvr_max--;
			break;
		    case 1:
			yrasmult = (float) ypix / (float) (xvr_max - xvr_min);
			xrasmult = (float) xpix / (float) (yvr_max - yvr_min);
			xvr_max--;
			yvr_min++;
			break;
		    case 2:
			xrasmult = (float) xpix / (float) (xvr_max - xvr_min);
			yrasmult = (float) ypix / (float) (yvr_max - yvr_min);
			xvr_min++;
			yvr_min++;
			break;
		    case 3:
		    default:
			yrasmult = (float) ypix / (float) (xvr_max - xvr_min);
			xrasmult = (float) xpix / (float) (yvr_max - yvr_min);
			xvr_min++;
			yvr_max--;
			break;
		}

		if (wantras && dev.smart_raster)
		{
		    rasterline  = sf_ucharalloc(xpix + 7 + 2);
		    rasterline2 = sf_ucharalloc(xpix + 7 + 2);
		    outraster = sf_ucharalloc2(xpix,ypix);

/*
 * See whether we want to dither or not.
 * Dither if monochrome device and dithering has been asked for.
 * It's up to the device to decide whether to actually do this.
 */
		    dither_it = dither && mono;

/*
 * If the device is color, but the raster happens to all be shades
 * of grey, let the device know so it can take advantage of the
 * fact. This is available to the routine via the external variable
 * "ras_allgrey".	joe & david 10/03/94
 */

		    ras_allgrey = YES;
		    if (!mono)
		    {
			if (c == VP_BIT_RASTER)
			{
			    for (ii=0; ii<2; ii++)
			    {
				if (
				    (color_set [ras_offset * ii][_RED] !=
				     color_set [ras_offset * ii][_GREEN]) ||
				    (color_set [ras_offset * ii][_GREEN] !=
				     color_set [ras_offset * ii][_BLUE])
				    )
				{
				    ras_allgrey = NO;
				    break;
				}
			    }
			}
			else
			{
			    for (ii=0; ii<256; ii++)
			    {
				if (
				    (color_set [ras_offset + ii][_RED] !=
				     color_set [ras_offset + ii][_GREEN]) ||
				    (color_set [ras_offset + ii][_GREEN] !=
				     color_set [ras_offset + ii][_BLUE])
				    )
				{
				    ras_allgrey = NO;
				    break;
				}
			    }
			}
		    }

		    /*
		     * Read in the Raster data for "Smart" devices, ie, those
		     * which can stretch (and dither) their own raster.
		     */
		    num_rep = 0;
		    for (yrast = 0; yrast < ypix; yrast++)
		    {
			/*
			 * Read in the next raster line, if we have a new one
			 */
			if (num_rep <= 0)
			{
			    num_rep = geth (pltin);
			    if (num_rep <= 0)
				ERR (FATAL, name, "Bad Raster line multiplier");

			    pos = 0;
			new_pat:num_pat = geth (pltin);
			    num_byte = geth (pltin);
			    if (num_pat <= 0 || num_byte <= 0 ||
				pos + num_pat * num_byte > xpix)
				ERR (FATAL, name, "Raster line not length promised");

			    if (num_pat > 1)
			    {
				if (dither_it)
				{
				    READ_RASTER (
					rasterline2[j] = GREY_MAP (ras_offset + (int) fgetc (pltin)),
					rasterline2[j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
					);
				}
				else
				{
				    READ_RASTER (
					rasterline2[j] = COLOR_MAP (ras_offset + (int) fgetc (pltin)),
					rasterline2[j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
					);
				}
				for (j = 0; j < num_pat; j++)
				{
				    for (k = 0; k < num_byte; k++)
				    {
					rasterline[pos] = rasterline2[k];
					pos++;
				    }
				}
			    }
			    else
			    {
				if (dither_it)
				{
				    READ_RASTER (
					rasterline[pos + j] = GREY_MAP (ras_offset + (int) fgetc (pltin)),
					rasterline[pos + j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
					);
				}
				else
				{
				    READ_RASTER (
					rasterline[pos + j] = COLOR_MAP (ras_offset + (int) fgetc (pltin)),
					rasterline[pos + j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
					);
				}
				pos += num_byte;
			    }

			    if (pos < xpix)
				goto new_pat;
			    if (pos != xpix)
				ERR (FATAL, name, "Raster line not length promised");
      
			}
			num_rep--;

			for (ii = 0; ii < xpix; ii++)
			{
			    outraster[yrast][ii] = rasterline[ii];
			}
		    }

		    /* Smart form */
		    dev.raster (xpix, ypix, xvr_min, yvr_min, xvr_max, yvr_max,
				outraster, ras_orient, dither_it);

		    free (rasterline);
		    free (rasterline2);
		    free (*outraster);
		    free (outraster);
		}
		else
		{
		    wlimit (xwmin, xwmax, &xvr_min, &xvr_max);
		    wlimit (ywmin, ywmax, &yvr_min, &yvr_max);

		    switch (ras_orient)
		    {
			case 0:
			    xvr_max++;
			    yvr_max++;
			    xr_max = xvr_max - xvru_min;
			    xr_min = xvr_min - xvru_min;
			    yr_max = yvr_max - yvru_min;
			    yr_min = yvr_min - yvru_min;
			    yru_max = yvru_max - yvru_min;
			    break;
			case 1:
			    xvr_max++;
			    yvr_min--;
			    xr_max = yvru_max - yvr_min;
			    xr_min = yvru_max - yvr_max;
			    yr_max = xvr_max - xvru_min;
			    yr_min = xvr_min - xvru_min;
			    yru_max = xvru_max - xvru_min;
			    break;
			case 2:
			    xvr_min--;
			    yvr_min--;
			    xr_max = xvru_max - xvr_min;
			    xr_min = xvru_max - xvr_max;
			    yr_max = yvru_max - yvr_min;
			    yr_min = yvru_max - yvr_max;
			    yru_max = yvru_max - yvru_min;
			    break;
			case 3:
			default:
			    xvr_min--;
			    yvr_max++;
			    xr_max = yvr_max - yvru_min;
			    xr_min = yvr_min - yvru_min;
			    yr_max = xvru_max - xvr_min;
			    yr_min = xvru_max - xvr_max;
			    yru_max = xvru_max - xvru_min;
			    break;
		    }
		    xru_min = 0;

		    if (yr_max < yr_min || xr_max < xr_min || !wantras)
		    {
			/*
			 * We need to read through all the raster stuff, even if
			 * we never use it.
			 */
			yr_max = yr_min;
			xr_max = xr_min;
		    }

		    rasterline  = sf_ucharalloc(xpix + 7 + 2);
		    rasterline2 = sf_ucharalloc(xpix + 7 + 2);

/*
 * See whether we need to dither or not.
 * Dither if monochrome device and dithering has been asked for.
 */
		    dither_it = dither && mono;

		    if (xr_max > xr_min)
		    {
			outraster2 = sf_ucharalloc2(xr_max - xr_min,1);

			if (dither_it)
			{
			    outraster = sf_ucharalloc2(xr_max - xr_min,1);
			}
			else
			{
			    outraster = outraster2;
			}
		    }
		    else
		    {
			outraster2 = NULL;
			outraster = NULL;
		    }

		    /*
		     * Read in the Raster data
		     */
		    lastrast = -1;
		    num_rep = 0;
		    for (i = yr_max - 1; i >= yr_min - 1; i--)
		    {
			yrast = (yru_max - 1 - i) * yrasmult;
			if (i == yr_min - 1)
			{
			    /*
			     * Assure that the last bit of unused raster, if any,
			     * is read. This last time through the loop is a
			     * "dummy".
			     */
			    yrast = ypix - 1;
			}

			for (ii = 0; ii < (yrast - lastrast); ii++)
			{
			    /*
			     * Read in the next raster line, if we have a new one
			     */
			    if (num_rep <= 0)
			    {
				num_rep = geth (pltin);
				if (num_rep <= 0){
				    ERR (FATAL, name, "Bad Raster line multiplier(2)");
				}
				pos = 0;
			    new_pat2:num_pat = geth (pltin);
				num_byte = geth (pltin);
				if (num_pat <= 0 || num_byte <= 0 ||
				    pos + num_pat * num_byte > xpix){
				    ERR (FATAL, name, "Raster line not length promised");
				}

				if ((ii + num_rep >= yrast - lastrast
				     && xr_max > xr_min && i >= yr_min) )
				{
				    if (num_pat > 1)
				    {
					if (dither_it)
					{
					    READ_RASTER (
						rasterline2[j] = GREY_MAP (ras_offset + (int) fgetc (pltin)),
						rasterline2[j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
						);
					}
					else
					{
					    READ_RASTER (
						rasterline2[j] = COLOR_MAP (ras_offset + (int) fgetc (pltin)),
						rasterline2[j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
						);
					}
					for (j = 0; j < num_pat; j++)
					{
					    for (k = 0; k < num_byte; k++)
					    {
						rasterline[pos] = rasterline2[k];
						pos++;
					    }
					}
				    }
				    else
				    {
					if (dither_it)
					{
					    READ_RASTER (
						rasterline[pos + j] = GREY_MAP (ras_offset + (int) fgetc (pltin)),
						rasterline[pos + j + jj] = GREY_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
						);
					}
					else
					{
					    READ_RASTER (
						rasterline[pos + j] = COLOR_MAP (ras_offset + (int) fgetc (pltin)),
						rasterline[pos + j + jj] = COLOR_MAP (ras_offset * ((ibyte & (001 << (7 - j))) != 0))
						);
					}
					pos += num_byte;
				    }
				}
				else
				{
/*
 * We're just going to turn right around and read another one,
 * So throw this away!
 */
				    if (c == VP_BYTE_RASTER)
				    {
					for (j = 0; j < num_byte; j++)
					{
					    fgetc (pltin);
					}
				    }
				    else
				    {
					for (j = 0; j < num_byte; j += 8)
					{
					    fgetc (pltin);
					}
				    }
				    pos += num_pat * num_byte;
				}
				if (pos < xpix)
				    goto new_pat2;
				if (pos != xpix){
				    ERR (FATAL, name, "Raster line not length promised");
				}

				if (wantras && xr_max > xr_min && i >= yr_min)
				{
				    for (j = xr_min; j < xr_max; j++)
				    {
					outraster2[0][j - xr_min] =
					    rasterline[(int) ((j - xru_min) * xrasmult)];
				    }
				}
			    }
			    num_rep--;
			}
			lastrast = yrast;

			if (xr_max > xr_min && i >= yr_min)
			{
			    if (dither_it)
			    {
				dithline (outraster2[0], outraster[0], 
					  xr_max - xr_min, 
					  yr_max - 1 - i, dither);
			    }
			    /* Dumb forms */
      
			    switch (ras_orient)
			    {
				case 0:
				    dev.raster (yr_max - 1 - i, yr_max - yr_min,
						xvr_min, i + yvru_min,
						xr_max - xr_min, ras_orient, outraster,
						0, 0);
				    break;
				case 1:
				    dev.raster (yr_max - 1 - i, yr_max - yr_min,
						xvru_min + i, yvr_max,
						xr_max - xr_min, ras_orient, outraster,
						0, 0);
				    break;
				case 2:
				    dev.raster (yr_max - 1 - i, yr_max - yr_min,
						xvr_max, yvru_max - i,
						xr_max - xr_min, ras_orient, outraster,
						0, 0);
				    break;
				case 3:
				    dev.raster (yr_max - 1 - i, yr_max - yr_min,
						xvru_max - i, yvr_min,
						xr_max - xr_min, ras_orient, outraster,
						0, 0);
				    break;
			    }
			}
		    }
		    free (rasterline);
		    free (rasterline2);
		    if (outraster2 != NULL)
		    {
			free (*outraster2);
			free (outraster2);
			if (dither_it) {
			    free (*outraster);
			    free (outraster);
			}
		    }
		    if (!wantras)
		    {
			xxx[0] = xvru_min;
			yyy[0] = yvru_min;
			xxx[1] = xvru_min;
			yyy[1] = yvru_max;
			xxx[2] = xvru_max;
			yyy[2] = yvru_max;
			xxx[3] = xvru_max;
			yyy[3] = yvru_min;
			drawpolygon (4, xxx, yyy);
		    }
		}
		break;
	    case VP_MESSAGE:	/* Text message */
		getvpstring ();
		message (MESG_READY,NULL);
		message (MESG_MESSAGE,NULL);
		message (MESG_TEXT, txbuffer);
		message (MESG_TEXT, CRLF);
		message (MESG_DONE,NULL);
		break;
	    case VP_SETDASH:
		npts = geth (pltin);
		if (npts > MAXDASH)
		{
		    ERR (FATAL, name, "Too complicated a dash line pattern.");
		}
		dashon = npts;
		k = 0;
		dashsum = 0.;
		for (ii = 0; ii < npts * 2; ii++)
		{
		    dashes[ii] = dashscale * (float) geth (pltin) / RPERIN;
		    if (dashes[ii] < 0.)
			ERR (FATAL, name, "Negative dash distance.");

		    if (dashes[ii] != 0.)
			k = 1;

		    dashsum += dashes[ii];
		}
		if (!k)
		    dashon = NO;
		dev.attributes (NEW_DASH, dashon, 0, 0, 0);
		break;
	    default:		/* error */
		ERR (FATAL, name,
		     "invalid VPLOT command decimal %d character %c",
		     (int) c, (char) c);
		break;
	}
    }
End_of_file:
    end_of_file();
    return;
}

void end_of_file(void)
{

/*
 * End of main while loop. Either fall out here or jump here when you hit
 * the end of the file somewhere.
 */

    dev.close (CLOSE_FLUSH);	/* force last vector out */

/*
 * End the group for this frame
 */
    group_number--;
    if (group_number != 0)
    {
	ERR (WARN, name,
	     "group left unclosed at end of file");
	group_number = 0;
    }
    else
    {
	if ((erase & FORCE_INITIAL) && (erase & DO_LITERALS))
	    dev.attributes (END_GROUP, group_number, EOF_ERASE_EGROUP, 0, 0);
	else
	    dev.attributes (END_GROUP, group_number, EOF_IGNORE_EGROUP, 0, 0);
    }

/*
 * Exit of dovplot
 */
    return;


}

/*
 * reset variables that can be affected by vplot commands when
 * processing multiple plots, and don't stay set across pages
 */
static void reset (void)
{
    int             ii, jj, kk;

    xwmin = xWmin;		/* plot window parameters defaulted */
    xwmax = xWmax;		/* to maximum size	 */
    ywmin = yWmin;
    ywmax = yWmax;

    reset_windows ();

    fat = fatmult * (float) (fatbase);
    dev.attributes (NEW_FAT, fat, 0, 0, 0);

    if (cur_color != COLOR_MAP (DEFAULT_COLOR))
    {
	need_devcolor = YES;
	cur_color = COLOR_MAP (DEFAULT_COLOR);
    }

    txalign.hor = TH_NORMAL;
    txalign.ver = TV_NORMAL;
    dev.attributes (NEW_ALIGN, txalign.hor, txalign.ver, 0, 0);

    ii = -1;
    jj = -1;
    kk = -1;
    if (dev.txfont != default_txfont)
    {
	dev.txfont = default_txfont;
	ii = dev.txfont;
    }
    if (dev.txprec != default_txprec)
    {
	dev.txprec = default_txprec;
	jj = dev.txprec;
    }
    if (dev.txovly != default_txovly)
    {
	dev.txovly = default_txovly;
	kk = dev.txovly;
    }
    dev.attributes (NEW_FONT, ii, jj, kk, 0);

    dashon = NO;
    dev.attributes (NEW_DASH, dashon, 0, 0, 0);

    overlay = default_overlay;
    dev.attributes (NEW_OVERLAY, overlay, 0, 0, 0);
}

void reset_windows (void)
/*< reset windows >*/
{
    extern int      xwmax_last, ywmax_last, xwmin_last, ywmin_last;

    if (xwmax != xwmax_last || ywmax != ywmax_last
	|| xwmin != xwmin_last || ywmin != ywmin_last)
	dev.attributes (SET_WINDOW, xwmin, ywmin, xwmax, ywmax);

    xwmin_last = xwmin;
    ywmin_last = ywmin;
    xwmax_last = xwmax;
    ywmax_last = ywmax;
}


static void outline_window (void)
{
    if (need_devcolor == YES || cur_color != COLOR_MAP (DEFAULT_COLOR))
    {
	dev.attributes (SET_COLOR, COLOR_MAP (DEFAULT_COLOR), 0, 0, 0);
	need_devcolor = NO;
    }
    dev.vector (xwmin, ywmin, xwmax, ywmin, 0, 0);
    dev.vector (xwmax, ywmin, xwmax, ywmax, 0, 0);
    dev.vector (xwmax, ywmax, xwmin, ywmax, 0, 0);
    dev.vector (xwmin, ywmax, xwmin, ywmin, 0, 0);
    if (cur_color != COLOR_MAP (DEFAULT_COLOR))
    {
	dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
	need_devcolor = NO;
    }
}

static void fill_background (void)
{
    int x[4], y[4];
    int cur_color_save;

    cur_color_save = cur_color;

    if (cur_color != COLOR_MAP (BACKGROUND_COLOR))
    {
        cur_color = COLOR_MAP (BACKGROUND_COLOR);
	need_devcolor = YES;
        /* drawpolygon will call update_color(); */
    }

    /* If the user requests a limited device plotting area, honor it. */
    x[0] = xWmin; y[0] = yWmin;
    x[1] = xWmin; y[1] = yWmax; 
    x[2] = xWmax; y[2] = yWmax;
    x[3] = xWmax, y[3] = yWmin;
    drawpolygon (4, x, y);
    
    if (cur_color != cur_color_save)
    {
        cur_color = cur_color_save;
	need_devcolor = YES;
    }
}

static void getvpstring (void)
{
    char           *txptr;
    int             ii;

    txptr = txbuffer;

    while (1)
    {
	ii = getc (pltin);

	if (ii == EOF)
	{
	    ERR (FATAL, name,
		 "Unexpected EOF encountered in Text string");
	}

	*txptr = ii;
	txptr++;

	if (ii == 0)
	{
	    break;
	}

	if ((txptr - txbuffer) == txbuflen)
	{
	    txbuffer = (char *) realloc (txbuffer, (unsigned) (txbuflen + TXBUFLEN));
	    txptr = txbuffer + txbuflen;
	    txbuflen += TXBUFLEN;
	}
    }
    return;
}

void drawpolygon (int npts, int *x, int *y)
/*< draw polygon >*/
{
    int             i, j;
    struct vertex  *vertex;
    static int      point;

    j = 0;
    if (npts > (vxbuflen - 1))
    {
	free ((char *) vxbuffer);
	vxbuffer =
	    (struct vertex *) sf_alloc (npts + 1,sizeof (struct vertex));
    }
    vertex = vxbuffer;
    xnew = x[j];
    ynew = y[j];
    j++;
    vertex->x = xnew;
    vertex->y = ynew;
    vertex->next = vertex + 1;
    vertex++;
    i = npts - 1;
    while (i--)
    {
	vertex->x = x[j];
	vertex->y = y[j];
	j++;
	if ((vertex->x != xnew) || (vertex->y != ynew))
	{
	    xnew = vertex->x;
	    ynew = vertex->y;
	    vertex->next = vertex + 1;
	    vertex->last = vertex - 1;
	    vertex++;
	    continue;
	}
	npts--;
    }
    vertex--;
    vxbuffer->last = vertex;
    vertex->next = vxbuffer;

    ipat = 0;
    point = cur_color;
    pat[ipat].patbits = &point;
    pat[ipat].xdim = 1;
    pat[ipat].ydim = 1;
    pat[ipat].xdim_orig = 1;
    pat[ipat].ydim_orig = 1;
    update_color ();
    if (npts > 2)
    {
	if (shade)
	    dev.area (npts, vxbuffer);
	else
	    vecoutline (vxbuffer);
    }
}

static int getpolygon (int npts)
{
    int             i;
    struct vertex  *vertex;

    if (npts > (vxbuflen - 1))
    {
	free ((char *) vxbuffer);
	vxbuffer =
	    (struct vertex *) sf_alloc ((unsigned) (npts + 1),
					sizeof (struct vertex));
    }
    vertex = vxbuffer;
    GETXY (xnew, ynew);
    vertex->x = xnew;
    vertex->y = ynew;
    vertex->next = vertex + 1;
    vertex++;
    i = npts - 1;
    while (i--)
    {
	GETXY (vertex->x, vertex->y);
	if ((vertex->x != xnew) || (vertex->y != ynew))
	{
	    xnew = vertex->x;
	    ynew = vertex->y;
	    vertex->next = vertex + 1;
	    vertex->last = vertex - 1;
	    vertex++;
	    continue;
	}
	npts--;
    }
    vertex--;
    vxbuffer->last = vertex;
    vertex->next = vxbuffer;
    return (npts);
}

static void update_color (void)
{
    if (need_devcolor)
    {
	dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
	need_devcolor = NO;
    }
}

static void add_a_cor (char *filename, int xcor, int ycor)
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

void getapoint (void)
/*< get a point >*/
{
    while (1)
    {
	xret = dev.xmax + 1;
	yret = dev.ymax + 1;
	if ((int) dev.getpoint (&xret, &yret) ||
	    (xret > dev.xmax - 5 && yret > dev.ymax - 5))
	    break;
	devtovpxy (xret, yret, &xret, &yret);
	add_a_cor (interact, xret, yret);
    }
}


void set_table(void){
    int             col_tab_no, red, green, blue, grey, dist, min_dist, best_col=0;
    float           red_f, green_f, blue_f, red_fm, green_fm, blue_fm;
    int             red_mask, green_mask, blue_mask;
    int             red_0, green_0, blue_0;
    int             red_7, green_7, blue_7;
    int             c,ii,k,i;
    /*
     * The logic here is a bit tricky. Basically, it goes like this:
     * If the device actually has a settable color of that number,
     * then reset the device's color as asked. Otherwise, take the
     * closest color that we have and map to that. Color 0 is the
     * background color, and is special. Merely "close" colors are
     * not mapped to Color 0... it has to be exact. Even devices with
     * NO settable colors are still sent colors 0 through 7, and
     * these colors are then always left set to their original
     * defaults. In this way you can handle terminals like the GIGI
     * which have colors but not SETTABLE colors. Also remember that
     * the device itself will THEN have to permute colors 0 through 7
     * on top of all this to correspond with Vplot's definitions!
     */
    col_tab_no = geth (pltin);
    if (col_tab_no > MAX_COL || col_tab_no < 0)
    {
	ERR (FATAL, name, "Bad color table value %d (%d max)",
	     col_tab_no, MAX_COL);
    }
    red = geth (pltin);
    green = geth (pltin);
    blue = geth (pltin);
    if (red > MAX_GUN || green > MAX_GUN || blue > MAX_GUN
	|| red < 0 || green < 0 || blue < 0)
    {
	ERR (FATAL, name,
	     "Bad color level in color %d (%d,%d,%d), %d max",
	     col_tab_no, red, green, blue, MAX_GUN);
    }
    
/*
 * Fun and games with color values
 */
    red_f = (float) red / (float) MAX_GUN;
    green_f = (float) green / (float) MAX_GUN;
    blue_f = (float) blue / (float) MAX_GUN;
    
    red_fm =
	redmap[3] + redmap[0] * red_f + redmap[1] * green_f + redmap[2] * blue_f;
    if (red_fm < 0.)
	red_fm = 0.;
    if (red_fm > 1.)
	red_fm = 1.;
    red_fm = pow (red_fm, redpow);
    red = ROUND (red_fm * MAX_GUN);
    
    green_fm =
	greenmap[3] + greenmap[0] * red_f + greenmap[1] * green_f + greenmap[2] * blue_f;
    if (green_fm < 0.)
	green_fm = 0.;
    if (green_fm > 1.)
	green_fm = 1.;
    green_fm = pow (green_fm, greenpow);
    green = ROUND (green_fm * MAX_GUN);
    
    blue_fm =
	bluemap[3] + bluemap[0] * red_f + bluemap[1] * green_f + bluemap[2] * blue_f;
    if (blue_fm < 0.)
	blue_fm = 0.;
    if (blue_fm > 1.)
	blue_fm = 1.;
    blue_fm = pow (blue_fm, bluepow);
    blue = ROUND (blue_fm * MAX_GUN);
    
    /*
     * Do color mapping (colormask option) Leave color 0 alone.
     */
    if (col_tab_no != 0)
    {
	red_mask = red;
	green_mask = green;
	blue_mask = blue;
	
/* Background color */
	red_0 = color_set[0][_RED];
	green_0 = color_set[0][_GREEN];
	blue_0 = color_set[0][_BLUE];
	
/* Opposite of background color */
	red_7 = MAX_GUN - color_set[0][_RED];
	green_7 = MAX_GUN - color_set[0][_GREEN];
	blue_7 = MAX_GUN - color_set[0][_BLUE];
	
	if (red == red_7 && green == green_7 && blue == blue_7)
	{
	    if (!colormask[3])
	    {
		red_mask = red_0;
		green_mask = green_0;
		blue_mask = blue_0;
	    }
	}
	else
	{
	    if (colormask[4])
	    {
		if (!colormask[0])
		    red_mask = red_0;
		if (!colormask[1])
		    green_mask = green_0;
		if (!colormask[2])
		    blue_mask = blue_0;
	    }
	    else
	    {
		if (colormask[0])
		{
		    green_mask = red_mask;
		    blue_mask = red_mask;
		}
		if (colormask[1])
		{
		    red_mask = green_mask;
		    blue_mask = green_mask;
		}
		if (colormask[2])
		{
		    red_mask = blue_mask;
		    green_mask = blue_mask;
		}
	    }
	}
	
	red = red_mask;
	green = green_mask;
	blue = blue_mask;
    }
    
    /*
     * Process the new color table value
     */

    if (col_tab_no < dev.num_col)
    {
	/*
	 * A settable color.
	 */
	dev.attributes (SET_COLOR_TABLE, col_tab_no,
			red, green, blue);
	color_set[col_tab_no][STATUS] = SET;
    }
    else if (col_tab_no > 7)
    {
	color_set[col_tab_no][STATUS] = MAPPED;
    }
    else
    {
	if (col_tab_no > 0)
	{
	    color_set[col_tab_no][STATUS] = MAP_SET;
	    /*
	     * Means that we need to map these but they are still set
	     * to the original colors and can't be changed. Color 0
	     * is always background, and any other color exactly the
	     * same color as color 0 "wants to be", even if color 0
	     * actually can't be and hasn't been reset, is mapped to
	     * color zero if it can't be set.
	     */
	}
    }
    /*
     * Save the color that this color table number wants to be
     */
    color_set[col_tab_no][_RED] = red;
    color_set[col_tab_no][_GREEN] = green;
    color_set[col_tab_no][_BLUE] = blue;
    /*
     * grey level is "Black and White TV style" mapped from color,
     * with corrections added by the subroutine greycorr.
     */
    grey = (blue * 1 + red * 2 + green * 4 + 6) / 7;
    color_set[col_tab_no][_GREY] = greycorr (grey);
    /*
     * If the next command is also a set color table command, we can
     * postpone doing this and kill 2 (or more) birds with one stone.
     */
    c = getc (pltin);
    if (c == EOF) end_of_file();
    ungetc ((char) c, pltin);
    if (c != VP_SET_COLOR_TABLE)
    {
	if (mono)
	{
	    for (ii = 1; ii <= MAX_COL; ii++)
	    {
		if (color_set[ii][_RED] == color_set[0][_RED] &&
		    color_set[ii][_GREEN] == color_set[0][_GREEN] &&
		    color_set[ii][_BLUE] == color_set[0][_BLUE])
		{
		    color_set[ii][MAP] = 0;
		}
		else
		{
		    color_set[ii][MAP] = 7;
		}
	    }
	}
	else
	{
	    /*
	     * For all color table entries that aren't set, (because
	     * we ran out of colors, for example) find the best color
	     * that IS set and use that instead.
	     */
	    for (ii = dev.num_col; ii <= MAX_COL; ii++)
	    {
		if (color_set[ii][STATUS] & MAPPED)
		{
		    min_dist = MAX_GUN * MAX_GUN * 8;
		    for (i = num_col_8 - 1; i >= 0; i--)
		    {
			/*
			 * Colors 1 through 7 are guaranteed SET, So
			 * we always get an answer. Color zero is
			 * background and special. To map to it you
			 * have to hit its color exactly, and no
			 * other color matched exactly first.
			 */
			if (color_set[i][STATUS] & SET)
			{
			    if (color_set[i][STATUS] == SET)
			    {
				k = color_set[i][_RED] - color_set[ii][_RED];
				dist = 2 * k * k;
				k = color_set[i][_GREEN] - color_set[ii][_GREEN];
				dist += 4 * k * k;
				k = color_set[i][_BLUE] - color_set[ii][_BLUE];
				dist += k * k;
			    }
			    else
			    {
				k = MAX_GUN * ((i & 2) / 2) - color_set[ii][_RED];
				dist = 2 * k * k;
				k = MAX_GUN * ((i & 4) / 4) - color_set[ii][_GREEN];
				dist += 4 * k * k;
				k = MAX_GUN * ((i & 1) / 1) - color_set[ii][_BLUE];
				dist += k * k;
			    }
			    if (dist < min_dist && (i != 0 || dist == 0))
			    {
				min_dist = dist;
				best_col = i;
				if (dist == 0)
				{
				    /*
				     * Might as well look no further
				     */
				    break;
				}
			    }
			}
		    }
		    color_set[ii][MAP] = best_col;
		}
	    }
	}
    }
}
