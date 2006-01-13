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
 *  source file:   ./filters/pslib/psarea.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), March 25 1988
 *      Wrote psarea.
 * Steve Cole (SEP), March 30 1988
 *      Routine now checks for "all black" and "all white" patterns and
 *      does the fill for these using setgray rather than messing with the
 *      halftone screen.
 *      Also, when the halftone screen has been altered, the routine now
 *      saves that state and avoids having to do all the work again if the
 *      same pattern is used to fill another area.
 * Steve Cole (SEP), March 31 1988
 *      Corrected orientation of pattern. Corrected scaling bug in
 *      call to setpattern.
 * Steve Cole (SEP), October 19 1988
 *      For mono=n case added a gsave to go with the grestore.
 * Joe Dellinger (Amoco), April 25 1996
 *	Use relative moves when they save space.
 */

#define mask0 ((unsigned char) (((unsigned int) 1) << 7))
#include <stdio.h>
#include <string.h>
#include "ps.h"
#include "../include/vertex.h"
#include "../include/extern.h"
#include "../include/err.h"
#include "../include/params.h"
#include "../include/pat.h"


extern FILE    *pltout;
extern int      ipat;
extern struct pat pat[];

psarea (npts, verlist)
    int             npts;
    struct vertex  *verlist;
{
unsigned char   mask, outchar;
int             nx, ny, nxold, nyold, i, ix, iy;
int            *bptr;
struct vertex  *v;
static int      ever_called = 0;
static int      error_size = 0;
static int      error_square = 0;
int             all_black, all_white;
static int      last_pattern = 0;
extern float    ps_curcolor;

int             ps_oldx = 0, ps_oldy = 0;
char            stringd[80];
char            stringr[80];

    endpath ();

    if (mono)
    {
/*
 * check for "all black" or "all white" patterns. It is a waste of time
 * to do these using the halftone screen method used below, since a
 * simple "setgray" call will work just as well.
 */
	nx = pat[ipat].xdim;
	ny = pat[ipat].ydim;
	bptr = pat[ipat].patbits;
	all_black = 1;
	all_white = 1;
	for (iy = 0; iy < ny; iy++)
	{
	    for (ix = 0; ix < nx; ix++)
	    {
		if (bptr[iy * nx + ix] != 0)
		    all_white = 0;
		if (bptr[iy * nx + ix] == 0)
		    all_black = 0;
	    }
	}
/*
 * This section of this routine is borrowed from the PostScript Language
 * Cookbook put out by Adobe Systems Incorporated, program 15 "Filling an Area
 * with a Pattern". I have commented the code to indicate the purpose
 * of each section and how the sections interact, but the Cookbook includes
 * more detailed comments that anyone attempting to modify this code
 * would certainly find useful.
 *
 * setuserscreendict allocates space for the routine setuserscreen.
 * ever_called is used to ensure that we only do all this initialization
 * once.
 * we don't need to do this initialization if the pattern is all black
 * or all white, because we use setgray rather than setscreen for those.
 */
	if (all_black == 0 && all_white == 0)
	{
	    if (ever_called == 0)
	    {
		fprintf (pltout, "/setuserscreendict 22 dict def\n");
		fprintf (pltout, "setuserscreendict begin\n");
		fprintf (pltout, " /tempctm matrix def\n");
		fprintf (pltout, " /temprot matrix def\n");
		fprintf (pltout, " /tempscale matrix def\n");
		fprintf (pltout, " /concatprocs\n");
		fprintf (pltout, " {/proc2 exch cvlit def\n");
		fprintf (pltout, "  /proc1 exch cvlit def\n");
		fprintf (pltout, "  /newproc proc1 length proc2 length add\n");
		fprintf (pltout, "   array def\n");
		fprintf (pltout, "  newproc 0 proc1 putinterval\n");
		fprintf (pltout, "  newproc proc1 length proc2 putinterval\n");
		fprintf (pltout, "  newproc cvx\n");
		fprintf (pltout, " } def\n");
		fprintf (pltout, "/resmatrix matrix def\n");
		fprintf (pltout, "/findresolution\n");
		fprintf (pltout, " {72 0 resmatrix defaultmatrix dtransform\n");
		fprintf (pltout, "  /yres exch def /xres exch def\n");
		fprintf (pltout, "  xres dup mul yres dup mul add sqrt\n");
		fprintf (pltout, " } def\n");
		fprintf (pltout, "end\n");
/*
 * setuserscreen takes the desired halftone cell size and angle in the current
 * user space from the procedure "setpattern" (defined below) and converts
 * them to device space values. Then it calls the postscript procedure
 * "setscreen" to re-define the halftone screen.
 */
		fprintf (pltout, "/setuserscreen\n");
		fprintf (pltout, " {setuserscreendict begin\n");
		fprintf (pltout, "  /spotfunction exch def\n");
		fprintf (pltout, "  /screenangle exch def\n");
		fprintf (pltout, "  /cellsize exch def\n");
		fprintf (pltout, "  /m tempctm currentmatrix def\n");
		fprintf (pltout, "  /rm screenangle temprot rotate def\n");
		fprintf (pltout, "  /sm cellsize dup tempscale scale def\n");
		fprintf (pltout, "  sm rm m m concatmatrix m concatmatrix pop\n");
		fprintf (pltout, "  1 0 m dtransform /y1 exch def /x1 exch def\n");
		fprintf (pltout, "  /veclength x1 dup mul y1 dup mul add sqrt def\n");
		fprintf (pltout, "  /frequency findresolution veclength div def\n");
		fprintf (pltout, "  /newscreenangle y1 x1 atan def\n");
		fprintf (pltout, "  m 2 get m 1 get mul m 0 get m 3 get mul sub\n");
		fprintf (pltout, "   0 gt\n");
		fprintf (pltout, "  {{neg} /spotfunction load concatprocs\n");
		fprintf (pltout, "    /spotfunction exch def\n");
		fprintf (pltout, "  } if\n");
		fprintf (pltout, "  frequency newscreenangle /spotfunction load\n");
		fprintf (pltout, "   setscreen\n");
		fprintf (pltout, " end\n");
		fprintf (pltout, " } def\n");
/*
 * setpatterndict allocates space for the procedure setpattern.
 * bitpatternspotfunction is an entry in this dictionary that defines
 * the pattern. setpattern will load this and pass it to setuserscreen.
 */
		fprintf (pltout, "/setpatterndict 18 dict def\n");
		fprintf (pltout, "setpatterndict begin\n");
		fprintf (pltout, " /bitison\n");
		fprintf (pltout, "  {/ybit exch def /xbit exch def\n");
		fprintf (pltout, "   /bytevalue bstring ybit bwidth mul xbit 8 idiv\n");
		fprintf (pltout, "    add get def\n");
		fprintf (pltout, "   /mask 1 7 xbit 8 mod sub bitshift def\n");
		fprintf (pltout, "   bytevalue mask and 0 ne\n");
		fprintf (pltout, "  } def\n");
		fprintf (pltout, "end\n");
		fprintf (pltout, "/bitpatternspotfunction\n");
		fprintf (pltout, " {setpatterndict begin\n");
		fprintf (pltout, " /y exch def /x exch def\n");
		fprintf (pltout, " /xindex x 1 add 2 div bpside mul cvi def\n");
		fprintf (pltout, " /yindex y 1 add 2 div bpside mul cvi def\n");
		fprintf (pltout, " xindex yindex bitison\n");
		fprintf (pltout, "  {/onbits onbits 1 add def 1}\n");
		fprintf (pltout, "  {/offbits offbits 1 add def 0}\n");
		fprintf (pltout, "  ifelse\n");
		fprintf (pltout, " end\n");
		fprintf (pltout, "} def\n");
/*
 * setpattern sets up the halftone screen given the parameters passed
 * from showpattern and sends this definition to setuserscreen.
 */
		fprintf (pltout, "/setpattern\n");
		fprintf (pltout, " {setpatterndict begin\n");
		fprintf (pltout, "  /cellsz exch def\n");
		fprintf (pltout, "  /angle exch def\n");
		fprintf (pltout, "  /bwidth exch def\n");
		fprintf (pltout, "  /bpside exch def\n");
		fprintf (pltout, "  /bstring exch def\n");
		fprintf (pltout, "  /onbits 0 def /offbits 0 def\n");
		fprintf (pltout, "  cellsz angle /bitpatternspotfunction load\n");
		fprintf (pltout, "   setuserscreen\n");
		fprintf (pltout, "  {} settransfer\n");
		fprintf (pltout, "  offbits offbits onbits add div setgray\n");
		fprintf (pltout, " end\n");
		fprintf (pltout, "} def\n");
		ever_called = 1;
	    }
/*
 * The following code uses the definitions made above to set the
 * halftone screen to represent the user's requested pattern.
 * If the previous area fill was also done with this pattern, however,
 * the necessary transfer function has been saved in the postscript
 * procedure "vplottransfer", and we can just re-load it.
 */
	    if (last_pattern == ipat && ipat != 0)
	    {
		fprintf (pltout, "vplottransfer settransfer\n");
	    }
	    else
	    {
/*
 * Here's where the vplot work begins.
 * Define the fill pattern. It must be square and it must be a
 * multiple of 16 bits wide.
 */
		nxold = pat[ipat].xdim;
		nx = nxold;
		if (nx % 16 != 0)
		{
		    nx = ((nx + 15) / 16) * 16;
		    if (error_size == 0)
		    {
			ERR (WARN, name, "Pspen pattern must be a multiple of 16 bits wide, so partially replicating");
			error_size = 1;
		    }
		}
		nyold = pat[ipat].ydim;
		ny = nyold;
		if (ny % 16 != 0)
		{
		    ny = ((ny + 15) / 16) * 16;
		    if (error_size == 0)
		    {
			ERR (WARN, name, "Pspen pattern must be a multiple of 16 bits wide, so partially replicating");
			error_size = 1;
		    }
		}
		if (nx != ny)
		{
		    if (error_square == 0)
		    {
			ERR (WARN, name, "Pspen pattern must be square, so partially replicating");
			error_square = 1;
		    }
		    if (nx > ny)
			ny = nx;
		    if (ny > nx)
			nx = ny;
		}
		bptr = pat[ipat].patbits;
		fprintf (pltout, "/pattern <");
		mask = mask0;
		outchar = 0;
		for (iy = 0; iy < ny; iy++)
		{
		    for (ix = 0; ix < nx; ix++)
		    {
			if (bptr[(nyold - 1 - (iy % nyold)) * nxold + (ix % nxold)] != 0)
			{
			    outchar |= mask;
			}
			if (mask == 1)
			{
			    fprintf (pltout, "%02.2x", outchar);
			    mask = mask0;
			    outchar = 0;
			}
			else
			{
			    mask >>= 1;
			}
		    }
		}

		fprintf (pltout, "> def\n");
/*
 * Save the graphics state (because we're going to destroy the halftone
 * screen, etc.) Construct the polygon and fill it by calling setpattern.
 * Then restore the graphics state.
 * arguments to setpattern:
 * number of bits on a side of the pattern
 * number of bytes on a side of the pattern
 * size of the pattern in PostScript units.
 */
		fprintf (pltout, "gsave\n");
		fprintf (pltout, "pattern %d %d 0 %d setpattern\n", nx, nx / 8, nx);
/*
 * Save the transfer function (this is the modified halftone screen)
 * so that if another area is filled with the same pattern we can avoid
 * all the work of computing the transfer function.
 */
		fprintf (pltout, "/vplottransfer currenttransfer def\n");
	    }
	}
/*
 * For all-black or all-white patterns, we simply use setgray.
 */
	else
	{
	    fprintf (pltout, "gsave\n");
	    if (all_white == 1 && ps_curcolor != PS_WHITE)
		fprintf (pltout, "1 setgray\n");
	    if (all_black == 1 && ps_curcolor != PS_BLACK)
		fprintf (pltout, "0 setgray\n");
	}
	fprintf (pltout, "np\n");

	v = verlist;
	fprintf (pltout, "%d %d m\n", v->x, v->y);

	ps_oldx = v->x;
	ps_oldy = v->y;
	v++;

	for (i = 1; i < npts; i++)
	{
	    /*
	     * A bit brute force, but a sure way of using whichever is
	     * shorter: absolute or relative draws.
	     */
	    sprintf (stringd, "%d %d d\n", v->x, v->y);
	    sprintf (stringr, "%d %d r\n", v->x - ps_oldx, v->y - ps_oldy);

	    if (strlen (stringd) <= strlen (stringr))
		fprintf (pltout, "%s", stringd);
	    else
		fprintf (pltout, "%s", stringr);

	    ps_oldx = v->x;
	    ps_oldy = v->y;
	    v++;
	}
	fprintf (pltout, "cf\n");
	fprintf (pltout, "gr\n");
    }
    else
    {
	fprintf (pltout, "gsave\n");
	fprintf (pltout, "np\n");

	v = verlist;
	fprintf (pltout, "%d %d m\n", v->x, v->y);

	ps_oldx = v->x;
	ps_oldy = v->y;
	v++;

	for (i = 1; i < npts; i++)
	{
	    /*
	     * A bit brute force, but a sure way of using whichever is
	     * shorter: absolute or relative draws.
	     */
	    sprintf (stringd, "%d %d d\n", v->x, v->y);
	    sprintf (stringr, "%d %d r\n", v->x - ps_oldx, v->y - ps_oldy);

	    if (strlen (stringd) <= strlen (stringr))
		fprintf (pltout, "%s", stringd);
	    else
		fprintf (pltout, "%s", stringr);

	    ps_oldx = v->x;
	    ps_oldy = v->y;
	    v++;
	}
	fprintf (pltout, "cf\n");
	fprintf (pltout, "gr\n");
    }

}
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
 *  source file:   ./filters/pslib/psattr.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), March 27 1988
 *      Routine no longer modifies the current vplot color. It modifies
 *      postscript's current color definition only. This is the way
 *      it's supposed to be done.
 * Steve Cole (SEP), March 29 1988
 *      Clipping (the SET_WINDOW case) is now supported.
 * Steve Cole (SEP), March 31 1988
 *      Dashed lines (NEW_DASH) is now supported.
 * Steve Cole (SEP), June 20 1988
 *      Dashon was incorrectly used as a local variable here. Cleaned
 *      up the clipping code some.
 * Steve Cole (SEP), August 19 1988
 *      Added newpath after setting of clipping window. Otherwise when
 *      a stroke command came along, the clipping path would be drawn in.
 * Joe Dellinger (SEP), October 19 1988
 *	Added "grey" option.
 * Steve Cole (SEP), November 1 1991
 *      Replaced grey option with proper color support. Colors are
 *      shaded in grey (like the grey option) by default in PostScript
 *      on greyscale devices.
 * Joe Dellinger (SOEST), October 15 1992
 *	After setting a new clipping window the global plot parameters
 *	line width, dash pattern, and color were being lost. (Polygon
 *	fill pattern is apparently not lost.)
 * Joe Dellinger (SOEST), October 16 1992
 *	Postscript allows through everything STRICTLY INSIDE a clipping
 *	path; the path itself is not included. Thus we have to extend the
 *	postscript clipping box by one minimal coordinate unit in each
 *	direction, because vplot assumes the edges of the clipping window
 *	ARE included.
 * Steve Cole (SEP), November 30 1992
 *      Added color raster support.
 * Joe Dellinger (Amoco), April 24 1996
 *	Wrap the stuff that needs to be done after a grestore up into
 *	a subroutine so it can be called from other routines too.
 * Joe Dellinger (Amoco), March 20 1997
 *	Don't do a grestore without a preceding gsave to go with it.
 * Joe Dellinger (BP Amoco), Oct 5 1999
 *	pspen color=y force=y should just give you the colors you ask for,
 *	period. pspen color=y force=n should give you their complements,
 *	period. Don't try to get excessively tricky flipping some colors
 *	around and not others. Some of the color versus monochrome logic
 *	was entangled.
 */

#include <stdio.h>
#include "./../include/attrcom.h"
#include "./../include/params.h"
#include "./../include/extern.h"
#include "./../include/enum.h"
#include "ps.h"

extern FILE    *pltout;
/* Is a dash pattern currently in effect? */
int             ps_dash_pattern_set = NO;
/* Is there a dash pattern saved in postscript's memory? */
int             ps_dash_pattern_exists = NO;
float           psscale;

/* black, blue, red, magenta, green, cyan, yellow, white, ... */
int             red[32768] = {0, 0, 255, 255, 0, 0, 255, 255};
int             green[32768] = {0, 0, 0, 0, 255, 255, 255, 255};
int             blue[32768] = {0, 255, 0, 255, 0, 255, 0, 255};

int             ps_grey_ras[32768];

float           ps_curcolor = 0.;	/* Only used for mono */
static float    new_ps_curcolor;	/* Only used for mono */
static int      ps_curcolor_no = 7;
static int      new_ps_curcolor_no;
static int      ps_curcolor_set = NO;

int             ps_done_clipping_gsave = NO;

psattributes (command, value, v1, v2, v3)
    register int    command, value;
    int             v1, v2, v3;
{
int             xmin, ymin, xmax, ymax;

/*
 * Are we in the middle of a move - draw - draw - draw ... sequence?
 * If so, wrap it up into a path and draw it, so we can move on and do
 * whatever needs to be done here!
 */
    endpath ();

    switch (command)
    {
    case SET_COLOR:
	ps_set_color (value);
	break;

    case SET_COLOR_TABLE:
/*
 * Note this will never be called if monochrome.
 */
	red[value] = v1;
	green[value] = v2;
	blue[value] = v3;

#ifdef PSDEBUG
	fprintf (stderr, "Color %d is now %d %d %d\n", value, v1, v2, v3);
#endif

/*
 * Used for greyscale raster
 */
	ps_grey_map (value);

/*
 * Check to see if the color they've just redefined is the one we're
 * currently using!
 */
	if (ps_curcolor_no == value)
	{
	    ps_curcolor_set = NO;
	    ps_set_color (ps_curcolor_no);
	}
	break;

/*
 * It is important that the clipping code begins with a grestore
 * and does a gsave just before setting the clipping path.
 * In PostScript, the clipping path is defined relative to the
 * current clipping path. So if you clipped with some path, then
 * tried to clip with a path that didn't overlap with the original
 * clipping region, you would get nothing. So we gsave before setting
 * the clipping window, then grestore just before re-setting it.
 * The first time through, there is no saved graphics state to restore.
 * Note: we have to check for an "inside out" clipping window
 * (xmin > xmax or ymin > ymax) here, because PostScript will not
 * catch this case.
 */
    case SET_WINDOW:
/*
 * We have to bump the box out by one unit in each direction, because
 * postscript doesn't include the boundary, but vplot does. At least, so says
 * the postscript reference manual. HOWEVER, the SPARC version of postscript
 * doesn't seem to actually do that: it includes the boundary at the
 * upper end, but doesn't at the lower end!!!! Under SPARC's interpretation
 * the clipping windows vary depending on the ORIENTATION of the same
 * postscript plot on the physical page, so I'm going to do it "the right way"
 * here. This allows an extra pixel width of stuff on each axis to sneak
 * through when plotting on the SPARC, unfortunately, but this seems the
 * best compromise to me. -- Joe D.
 */
	xmin = value - 1;
	ymin = v1 - 1;
	xmax = v2 + 1;
	ymax = v3 + 1;

/*
 * Inside-out clipping window test. We use (value,v1) and (v2,v3) instead
 * of (xmin,ymin) and (xmax,ymax), because we need to do the test according
 * to vplot's idea of how clipping windows work, not Postscript's.
 */
	if (value > v2)
	    xmin = xmax;
	if (v1 > v3)
	    ymin = ymax;

	if (ps_done_clipping_gsave == YES)
	    fprintf (pltout, "grestore\n");

	fprintf (pltout, "/clippingpath\n");
	fprintf (pltout, " {newpath\n");
	fprintf (pltout, " %d %d m\n", xmin, ymin);
	fprintf (pltout, " %d %d d\n", xmax, ymin);
	fprintf (pltout, " %d %d d\n", xmax, ymax);
	fprintf (pltout, " %d %d d\n", xmin, ymax);
	fprintf (pltout, " closepath} def\n");

	fprintf (pltout, "gsave clippingpath eoclip\n");
	ps_done_clipping_gsave = YES;

	fprintf (pltout, "newpath\n");

	/*
	 * The grestore caused postscript to forget most of the current
	 * state. Set it all back up again.
	 */
	ps_fixup_after_grestore ();

	/* Global plot parameters fixed; we can exit now. */
	break;

    case NEW_DASH:
	ps_set_dash (value);
	break;

    default:
	break;
    }

    return 0;
}


ps_fixup_after_grestore ()
{
/*
 * Alas, when we did the "gerestore" above ALL GLOBAL PLOT SETTINGS WERE LOST!
 * So go back and tidy up. The current fatness is easy; we just need to let
 * psvector know that the current value is the default (0) again. As for the
 * dash pattern, the variable ps_dash_pattern_exists tells psvector it's been
 * lost and the routine ps_set_dash provides a way for psvector to recreate it
 * when needed. The current color has to be reset right away.
 */

    /* The current fatness reverts to 0 */
    ps_last_fat = 0;

    /* No dash pattern is in effect */
    ps_dash_pattern_set = NO;
    /* No dash pattern remains in memory */
    ps_dash_pattern_exists = NO;

    /*
     * The current color got reset to the default. We could try to be tricky
     * and not reset the color if the default is OK, but it's safer and
     * simpler to force the color to be explicitly reset right away.
     */
    ps_curcolor_set = NO;
    ps_set_color (ps_curcolor_no);
}

ps_grey_map (coltab)
    int             coltab;
{
int             grey;

    /*
     * Calculate the grey level that goes with this color.
     */
    grey = (blue[coltab] * 1 + red[coltab] * 2 + green[coltab] * 4 + 6) / 7;

    /*
     * For raster, apply appropriate gamma corrections and invert as
     * necessary.
     */
    ps_grey_ras[coltab] = greycorr (grey);
}

ps_set_color (value)
    int             value;
{
    new_ps_curcolor_no = value;


#ifdef PSDEBUG
    fprintf (stderr, "Set? %d, Change Color to %d from %d\n",
	     ps_curcolor_set, new_ps_curcolor_no, ps_curcolor_no);
#endif


    if (mono)
    {
	if (new_ps_curcolor_no > 0)
	    new_ps_curcolor = PS_BLACK;
	else
	    new_ps_curcolor = PS_WHITE;

	if (new_ps_curcolor != ps_curcolor || !ps_curcolor_set)
	{
	    fprintf (pltout, "%.2g setgray\n", new_ps_curcolor);

	    ps_curcolor = new_ps_curcolor;
	    ps_curcolor_no = new_ps_curcolor_no;
	    ps_curcolor_set = YES;
	}
    }
    else
    {
	if (new_ps_curcolor_no != ps_curcolor_no || !ps_curcolor_set)
	{
	    if (force_color)
	    {
#ifdef PSDEBUG
		fprintf (stderr, "%.2g %.2g %.2g setrgbcolor\n", red[new_ps_curcolor_no] / 255., green[new_ps_curcolor_no] / 255., blue[new_ps_curcolor_no] / 255.);
#endif
		fprintf (pltout, "%.2g %.2g %.2g setrgbcolor\n", red[new_ps_curcolor_no] / 255., green[new_ps_curcolor_no] / 255., blue[new_ps_curcolor_no] / 255.);
	    }
	    else
	    {
#ifdef PSDEBUG
		fprintf (stderr, "%.2g %.2g %.2g setrgbcolor\n", 1. - red[new_ps_curcolor_no] / 255., 1. - green[new_ps_curcolor_no] / 255., 1. - blue[new_ps_curcolor_no] / 255.);
#endif
		fprintf (pltout, "%.2g %.2g %.2g setrgbcolor\n", 1. - red[new_ps_curcolor_no] / 255., 1. - green[new_ps_curcolor_no] / 255., 1. - blue[new_ps_curcolor_no] / 255.);
	    }

	    ps_curcolor_no = new_ps_curcolor_no;
	    ps_curcolor_set = YES;
	}
    }
}

ps_set_dash (value)
    int             value;
{
/*
 * A local variable named dashon. This routine has no business resetting
 * the global variable of the same name.
 */
int             dashon;
int             ii;

    dashon = value;

/*
 * No need to erase the pattern; psvector knows not to use it if dashon==0.
 */
    if (dashon == 0)
	return;

    fprintf (pltout, "[");
    for (ii = 0; ii < dashon * 2; ii++)
    {
	fprintf (pltout, "%.2f ", dashes[ii] * PSPERIN / psscale);
    }
    fprintf (pltout, "] 0 setdash\n");

    /* A dash pattern is now in postscript's memory. */
    ps_dash_pattern_exists = YES;
    /* And it is currently active, too. */
    ps_dash_pattern_set = YES;
}
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
 *  source file:   ./filters/pslib/psclose.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), March 29 1988
 *	Moved showpage to pserase.
 * Joe Dellinger (SOEST), Feb 18 1992
 *	Added option for "PSPRINTER" environment variable.
 * Steve Cole (SEP), October 22 1992
 *	Added color printer support.
 * Joe Dellinger (AMOCO), March 22 1995
 *      Added custom paper size support.
 *      Added comment for Amoco users so they can see where their plot
 *      is being spooled. If you don't like that, #define NOCOMMENTS
 * Hector Urdaneta (SEP), April 26 1996
 *      Added ifdef SGI64, define NOCOMMENTS, endif 
 *      Not compiling in the SGI
 */

/* Routine to finish up */


#include <stdio.h>
#if defined(__stdc__) || defined(__STDC__)
#include <string.h>
#else
#include <strings.h>
#endif
#if defined(__stdc__) || defined(__STDC__)
#include <stdlib.h>
#else
extern char    *getenv ();
#endif

#include "../include/closestat.h"
#include "../include/enum.h"
#include "../include/extern.h"
#include "../include/err.h"
#include "ps.h"

#ifdef SEP
#define GETPAR getch
#else
#define GETPAR getpar
#endif

psclose (status)
    int             status;
{
char            system_call[120];
char            printer[40];
char           *stringptr;
int             ecode;

    endpath ();

    switch (status)
    {
    case CLOSE_NORMAL:
	/*
	 * if(tex == YES) fprintf(pltout, "%c", POP);
	 */

#ifndef NOCOMMENTS
	if ((strcmp (psprintertype, "default") != 0)
	    || (pixels_per_inch != DEFAULT_PIXELS_PER_INCH))
	{
	    ERR (COMMENT, name,
		 "Printer is of type \"%s\", with %g pixels per inch.",
		 psprintertype, pixels_per_inch);
	}
#endif

	/*
	 * If we created a temporary file, spool it.
	 */
	if (file_created)
	{
	    fclose (pltout);

	    if ((stringptr = getenv ("PSPRINTER")) != NULL)
		strcpy (printer, stringptr);
	    else if (mono)
		strcpy (printer, "postscript");
	    else
		strcpy (printer, "colorps");

	    GETPAR ("printer", "s", printer);

	    if (ps_set_papersize)
	    {
		sprintf (system_call,
			 "lpr -r -s -P%s -CX%.2f,Y%.2f %s",
			 printer, ps_ylength, ps_xlength, scratch_file);
	    }
	    else
	    {
		sprintf (system_call,
			 "lpr -r -s -P%s %s", printer, scratch_file);
	
	    }

#ifndef NOCOMMENTS
	    if (ps_set_papersize)
	    {
		ERR (COMMENT, name,
		     "Spooling plot using command \"%s\".", system_call);
	    }
	    else
	    {
		ERR (COMMENT, name,
		     "Spooling plot to printer \"%s\".", printer);
	    }
#endif
	    ecode = system (system_call);

	    /*
	     * Shift 8 bits over; what's left contains the exit code from lpr
	     */
	    if ((ecode & 0xFF) != 0)
	    {
		ERR (WARN, name,
		"Signal stopped or killed system call \"%s\".", system_call);
		ERR (WARN, name,
		     "Output postscript file may have been left behind in \"%s\".",
		     scratch_file);
	    }
	    else
	    {
		ecode = (ecode > 8);
		if (ecode != 0)
		{
		    ERR (WARN, name,
			 "Exit code %d from lpr. Is \"%s\" a valid printer?",
			 ecode, printer);
		    if (stringptr == NULL)
		    {
			ERR (COMMENT, name,
			     "Perhaps you need to do \"setenv PSPRINTER printer_name\"?");
		    }
		    else
		    {
			ERR (COMMENT, name,
			     "Is your environment variable $PSPRINTER = \"%s\" correct?",
			     printer);
		    }
		    ERR (WARN, name,
			 "The output postscript file may have been left behind in \"%s\".",
			 scratch_file);
		}
	    }
	}
	file_created = NO;
	break;
    case CLOSE_ERROR:
    case CLOSE_NOTHING:
    case CLOSE_INTERRUPT:
	break;
    case CLOSE_DONE:
	/*
	 * If we created a temporary file, remove it
	 */
	if (file_created)
	    unlink (scratch_file);
	break;
    case CLOSE_FLUSH:
	fflush (pltout);
	break;
    case CLOSE_PAUSE:
	break;
    default:			/* not meant for us, ignore */
	break;
    }
}
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
 *  source file:   ./filters/pslib/psconf.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), August 26 1987
 *	Replaced genraster with psraster.
 * Steve Cole (SEP), March 25 1988
 *      Added psarea as it is now available.
 */

/*
 * Apple Laser Writer Configuration
 */
#include <stdio.h>
#include "../include/enum.h"
#include "../include/extern.h"

/*
 * mandatory declarations and initializations
 */
#ifdef SEP
char            name[] = "Pspen";
#else
char            name[] = "pspen";
#endif
#include "psdoc.h"

/*
 * device routine table
 */
extern int
psopen (), psreset (), pserase (), psclose ();
extern int
psvector (), psplot (), genmarker (), pstext (), psarea ();
extern int
smart_psraster (), psattributes ();
extern int
genstartpoly (), genmidpoly (), genendpoly ();
extern int
nulldev (), genpoint (), genmessage (), genmarker ();

struct device   dev = {

 /* control routines */
		       psopen,	/* open */
		       psreset,	/* reset */
		       genmessage,	/* message */
		       pserase,	/* erase */
		       psclose,	/* close */

 /* high level output */
		       psvector,	/* vector */
		       genmarker,	/* marker */
		       pstext,	/* text */
		       psarea,	/* area */
		       smart_psraster,	/* raster */
		       genpoint,	/* point */
		       psattributes,	/* attributes */

 /* input */
		       nulldev,	/* getpoint */
		       nulldev,	/* interact */

 /* low level output */
		       psplot,	/* plot */
		       nulldev,	/* startpoly */
		       nulldev,	/* midpoly */
		       nulldev,	/* endpoly */
};
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
 *  source file:   ./filters/pslib/pserase.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), March 29 1988
 *      Added plot label.
 * Steve Cole (SEP), March 31 1988
 *      Added support for ncopies and hold command line options.
 * Dave Nichols (SEP), May 22 1990
 *	Added a grestore at the end of the page to match the gsave at the start
 *	Modified output for tex=y to be compatible with /psfig macros in TeX.
 * Joe Dellinger (SOEST), Oct 15 1992
 *	Added a grestore/gsave pair just before printing the label
 *	to turn off any clipping window which may be in effect.
 * Joe Dellinger (AMOCO), April 24 1996
 *	Paper size was not set on pages after the first.
 *	Line join type was not set if tex=y, causing "fanged S syndrome".
 *	The variable "rotate" should not have been used here.
 *	Fatness (and some other variables) needed to be reset at the
 *	beginning of each new page.
 * Joe Dellinger (AMOCO), March 20 1997
 *	If we did a gsave for a clipping window, then do a grestore
 *	here so gsaves and grestores are properly paired.
 */

#include	<stdio.h>
#include	"ps.h"
#include	"../include/erasecom.h"
#include	"../include/enum.h"
#include	"../include/extern.h"

/*
 * Location and height of label and surrounding box in ps coordinates
 */
#define TEXT_YPOS	100
#define TEXT_HEIGHT	50
#define TEXT_PAD	18

extern FILE    *pltout;
float           psscale;

extern int      ps_done_clipping_gsave;


pserase (command)
    int             command;
{
char            full_label[100];
static int      page_count = 1;

    endpath ();

    switch (command)
    {
    case ERASE_MIDDLE:
    case ERASE_END:
/*
 * put on label if desired
 */
	if (label[0] != '\0')
	{
	    if (page_count == 1 && command == ERASE_END)
	    {
		sprintf (full_label, "%s", label);
	    }
	    else
	    {
		sprintf (full_label, "%s : Page %d.", label, page_count);
	    }
	    /*
	     * Do a grestore so any clipping window currently in force won't
	     * clip away the plot label. Also do a gsave, because by
	     * convention the plot ends with a grestore.
	     */
	    fprintf (pltout, "grestore gsave\n");
	    fprintf (pltout, "/labelshow\n");
	    fprintf (pltout, " {dup stringwidth pop\n");
	    fprintf (pltout, "  dup dup %d exch sub\n", dev_xmax - 50);
	    fprintf (pltout, "  newpath %d m\n", TEXT_YPOS);
	    fprintf (pltout, "  gsave 0 %d rlineto\n", -TEXT_PAD);
	    fprintf (pltout, "  0 rlineto\n");
	    fprintf (pltout, "  0 %d rlineto\n", TEXT_HEIGHT + TEXT_PAD);
	    fprintf (pltout, "  neg 0 rlineto\n");
	    fprintf (pltout, "  closepath 1 setgray fill grestore\n");
	    fprintf (pltout, "  gsave 0 setgray\n");
	    fprintf (pltout, "  show grestore} def\n");
	    fprintf (pltout, "/Helvetica-BoldOblique findfont %d scalefont setfont\n", TEXT_HEIGHT);
	    fprintf (pltout, "(%s) labelshow\n", full_label);
	}
	if (tex == NO)
	{
	    if (ncopies_document != 1)
	    {
		fprintf (pltout, "/#copies %d def\n", ncopies_document);
	    }
	    if (hold == YES)
	    {
		fprintf (pltout, "statusdict begin\n");
		fprintf (pltout, "/manualfeed true def\n");
		fprintf (pltout, "/manualfeedtimeout 300 def\n");
		fprintf (pltout, "end\n");
	    }
	}

        if (ps_done_clipping_gsave == YES)
            fprintf (pltout, "grestore\n");
        ps_done_clipping_gsave = NO;

	fprintf (pltout, "grestore\n");
	fprintf (pltout, "showpage\n");

	if (command == ERASE_MIDDLE)
	{
	    if (tex == NO)
	    {
		fprintf (pltout, "initgraphics 1 setlinecap 1 setlinejoin\n");
		fprintf (pltout, "%d rotate", 90);
		fprintf (pltout, " 0 %.2f translate %.2f %.2f scale gsave\n",
			 -ps_ypapersize * PSPERIN, psscale, psscale);
	    }
	    else
	    {
		fprintf (pltout, " %.2f %.2f scale gsave\n", psscale, psscale);
		fprintf (pltout, " 1 setlinecap 1 setlinejoin\n");
	    }

	    /*
	     * Postscript has forgetten most of the current state. Set it all
	     * back up again.
	     */
	    ps_fixup_after_grestore ();
	}
	page_count++;
	break;
    case ERASE_START:
    default:
	break;
    }
}
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
 *  source file:   ./filters/pslib/psopen.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), September 15 1987
 *	Added smart_raster, mono, dither, pixc, and greyc to attributes.
 * Steve Cole (SEP), January 10 1988
 *      Vplot pen font is now default font.
 * Steve Cole (SEP), March 23 1988
 *      Default font is now vplot default hardcopy font.
 * Steve Cole (SEP), March 29 1988
 *      Added an intial gsave so that the first clipping window call
 *      will not destroy the rotation and translation done here when
 *      it tries to back out the previous clipping call.
 *	Added need_end_erase for plotting label.
 *      Corrected definition of limits of device space.
 * Steve Cole (SEP), March 31 1988
 *      The vplot to PostScript units conversion factor is computed here
 *	instead of being hardwired in many different places.
 *      Scaling is done here so that all coordinates output by pspen
 *      can be integers.
 * Joe Dellinger (SEP), October 19 1988
 *	Added "grey" and "harddither" options, fixed bugs, use "move/draw"
 *	method of sending paths.
 * Joe Dellinger (SEP), January 19 1989
 *	This is a "page type" hardcopy device;
 *	size=relative makes more sense.
 * Joe Dellinger (SEP), April 17 1990
 *	Remove pointless calloc call; "getlogin" provides its own allocation.
 * Dave Nichols (SEP), May 22 1990
 *	Changed output for tex=y to be compatible with /psfig macros in TeX.
 * Dave Nichols (SEP), Aug 14 1990
 *	Changed size of page to be maximum possible and force size=ABSOLUTE
 *	for tex=y to prevent clipping of plots off the page.
 * Steve Cole (SEP), Nov 1 1991
 *	Added color support.
 * Dave Nichols (SEP) May 2 1992
 *      Allow VPLOTSPOOLDIR environment variable to override PEN_SPOOL
 * Joe Dellinger (SOEST) Oct 15 1992
 *	The Apple LaserWriter page limits are a little more generous
 *	than the SPARC printer allows.
 * Steve Cole (SEP), Nov 30 1992
 *      Added color raster support.
 *	Added limits for Tektronix color postscript printer. These are
 *      used if pspen is invoked by the alias tcprpen, or if
 *      wstype=tcpr.
 * Ray Abma (SEP) 3 Aug 1993
 *      changed getdate to dateget to avoid conflict with HP time.h
 * David Lumley (SEP) March 7 1994
 *	added translate call so tcpr printable area doesn't get clipped
 * Joe Dellinger (AMOCO) March 22 1995
 *	Added the ability to specify a page size from the command line.
 *	Added the ability to set the default page size either by the
 *	define DEFAULT_PAPER_SIZE at compile time or by setting an
 *	environmental variable of the same name.
 *	Corrected the bug (I hope?) that caused plots to sometimes be
 *	lost if pspen was run in background from a shell after the user
 *	logged out. Added more detail to error messages.
 *	Added "oyo" option to wstype.
 * Joe Dellinger (AMOCO) April 24 1996
 *	Set line join type when tex=y too.
 *	The variable "rotate" was always zero because it hadn't been
 *	properly set yet and shouldn't have been used here.
 * Joe Dellinger (AMOCO) April 25 1996
 *	Add "s" and "r" as shorthands for "stroke" and "rlineto", respectively.
 * Joe Dellinger (AMOCO) April 26 1996
 *	Add "x" as shorthand for "0 0 rlineto" (a dot).
 * Joe Dellinger (BP Amoco) Oct 5, 1999
 *	Set the actual colors when calling psattributes (changes made
 *	there to fix a bug required a change here).
 */

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "../include/vplot.h"
#include "../include/err.h"
#include "../include/enum.h"
#include "../include/extern.h"
#include "../include/params.h"
#include "../include/attrcom.h"
#include "ps.h"
#include <pwd.h>

#if defined(__stdc__) || defined(__STDC__)
#include <stdlib.h>
#else
extern char    *getenv ();
#endif

#ifdef SEP
#define GETPAR getch
#else
#define GETPAR getpar
#endif


#ifndef DEFAULT_PAPER_SIZE
#define DEFAULT_PAPER_SIZE	"letter"
#endif

#ifndef DEFAULT_PAPER_UNITS
#define DEFAULT_PAPER_UNITS	'i'
#endif


int             file_created = NO;
int             tex = NO;
int             ncopies_document = 1;
char            label[100];
char            ps_paper[80];

char            mapfile[100] = "default";
char            scratch_file[100];
int             hold;

int             default_ps_font = DEFAULT_HARDCOPY_FONT;
float           psscale;

int             ps_color = NO;
int             dumb_fat = NO;
int             force_color = NO;
extern int      red[], green[], blue[];

extern int      smart_psraster (), dumb_psraster ();

float           ps_xlength, ps_ylength;
float           ps_xmin, ps_xmax, ps_ymin, ps_ymax;
float           ps_xborder, ps_yborder;
int             ps_set_papersize = NO;
float           ps_ypapersize;
char            psprintertype[80] = "default";

psopen ()
{
char           *user_name, *getlogin ();
char            date[50];
int             hard_dither = YES;
int             i;
char            units;
float           unitscale;
extern int      isafile ();
struct stat     statbuf;
int             creatafile;
int             yesget;
char           *paperstring;


/*
 * A list of different printer types recognized
 */
    if ((strcmp (callname, "oyopen") == 0) ||
	(strcmp (callname, "Oyopen") == 0) ||
	(strcmp (wstype, "oyo") == 0))
    {
	strcpy (psprintertype, "oyo");
    }
    else if ((strcmp (callname, "tcprpen") == 0) ||
	     (strcmp (callname, "Tcprpen") == 0) ||
	     (strcmp (wstype, "tcpr") == 0))
    {
	strcpy (psprintertype, "tcpr");
    }
    else
    {
/* Default */
	strcpy (psprintertype, "default");
    }


/*
 * physical device parameters
 */

/*
 * pixc and greyc are used even if smart_raster = YES because
 * PostScript's built in dithering does not compensate for the
 * effects that pixc and greyc correct for.
 */

    if ((strcmp (callname, "oyopen") == 0) ||
	(strcmp (callname, "Oyopen") == 0) ||
	(strcmp (wstype, "oyo") == 0))
    {
/*
 * OYO electrostatic plotters need a VERY different pixc and greyc
 */
	pixc = 1.5;
	greyc = .1;

/*
 * Our main OYO happens to be 300 dpi
 */
	pixels_per_inch = 300.;

/*
 * Change this if for example you want the oyo-s to use bigger sheets
 * of paper by default than the usual HP printer engines.
 */
	strcpy (ps_paper, DEFAULT_PAPER_SIZE);
    }
    else
    {
/* For HP printer engines, the usual default */

	pixc = 0.6;
	greyc = -0.5;

	pixels_per_inch = DEFAULT_PIXELS_PER_INCH;

	strcpy (ps_paper, DEFAULT_PAPER_SIZE);
    }


/*
 * For now, elsewhere in the code we just ignore aspect_ratio,
 * implicitly assuming it will always be 1.
 */
    aspect_ratio = 1.0;


/*
 * Allow override of pixels_per_inch from the command line. (If we used
 * GETPAR here it would also search the incoming history file... I think
 * it's safer to require them to be on the command line.)
 */
    getpar ("ppi", "f", &pixels_per_inch);


/*
 * based on the number of pixels per inch and PostScript's units
 * of length, we compute a PostScript scale factor here that will allow us to
 * send vplot coordinates (integers) directly to PostScript.
 */
    psscale = PSPERIN / ((float) pixels_per_inch);


/*
 * American paper sizes
 *
 * letter: 8.5 by 11 inch page; imageable region of 8.0 by 10.92 inches.
 * legal: 8.5 by 13 inch page; imageable region of 6.72 by 13 inches.
 * (these limits are imposed by the apple laserwriter).
 */


/*
 * Default units for custom paper sizes.
 */
    units = DEFAULT_PAPER_UNITS;

/*
 * Environment variable can override default paper size...
 */
    if ((paperstring = getenv ("DEFAULT_PAPER_SIZE")) != NULL)
	strcpy (ps_paper, paperstring);

/*
 * but GETPAR has the very last say.
 */
    if (GETPAR ("paper", "s", ps_paper))
	paperstring = NULL;

    if (strcmp (ps_paper, "letter") == 0)
    {
	if ((strcmp (callname, "tcprpen") == 0) ||
	    (strcmp (callname, "Tcprpen") == 0) ||
	    (strcmp (wstype, "tcpr") == 0))
	{
/*
 * Empirically determined plotting limits for the Tek color plotter with
 * letter-size paper
 */
	    dev_xmin = 1.23333 * pixels_per_inch + .5;
	    dev_ymin = .18333 * pixels_per_inch + .5;
	    dev_xmax = 9.76666 * pixels_per_inch + .5;
	    dev_ymax = 8.31333 * pixels_per_inch + .5;
/*
 * ps_ypapersize controls an overall plot shift. Normally it will be the
 * same as the height (Y) dimension of the paper, but may be different in
 * case the plotter doesn't put (0,0) at the corner of the paper.
 */
	    ps_ypapersize = 7.66666;
	}
	else
	{
/* Empirically determined limits for SPARC printers with letter-size paper */
	    dev_xmin = .19 * pixels_per_inch + .5;
	    dev_ymin = .26666 * pixels_per_inch + .5;
	    dev_xmax = 10.88333 * pixels_per_inch + .5;
	    dev_ymax = 8.23333 * pixels_per_inch + .5;
	    ps_ypapersize = 8.5;
	}
    }
    else if (strcmp (ps_paper, "legal") == 0)
    {
/* Empirically determined limits for SPARC printers with legal-size paper */
	dev_xmin = .50333 * pixels_per_inch + .5;
	dev_ymin = .89 * pixels_per_inch + .5;
	dev_xmax = 13.49666 * pixels_per_inch + .5;
	dev_ymax = 7.61 * pixels_per_inch + .5;
	ps_ypapersize = 8.5;
    }
    else if (strcmp (ps_paper, "a5") == 0 || strcmp (ps_paper, "A5") == 0)
    {
	/*
	 * For our European friends, A5 paper
	 * 
	 * For now, use (DEFAULT_BORDER * CMPERIN) as the best guess for a
	 * nonimageable area at the edges of the paper all around.
	 */

	dev_xmin = (DEFAULT_BORDER * CMPERIN) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_ymin = (DEFAULT_BORDER * CMPERIN) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_xmax = (21.0 - (DEFAULT_BORDER * CMPERIN)) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_ymax = (14.85 - (DEFAULT_BORDER * CMPERIN)) * (1. / CMPERIN) * pixels_per_inch + .5;
	ps_ypapersize = 14.85 * (1. / CMPERIN);
    }
    else if (strcmp (ps_paper, "a4") == 0 || strcmp (ps_paper, "A4") == 0)
    {
	/*
	 * For our European friends, A4 paper
	 * 
	 * For now, use (DEFAULT_BORDER * CMPERIN) as the best guess for a
	 * nonimageable area at the edges of the paper all around.
	 */

	dev_xmin = (DEFAULT_BORDER * CMPERIN) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_ymin = (DEFAULT_BORDER * CMPERIN) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_xmax = (29.7 - (DEFAULT_BORDER * CMPERIN)) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_ymax = (21.0 - (DEFAULT_BORDER * CMPERIN)) * (1. / CMPERIN) * pixels_per_inch + .5;
	ps_ypapersize = 21.0 * (1. / CMPERIN);
    }
    else if (strcmp (ps_paper, "a3") == 0 || strcmp (ps_paper, "A3") == 0)
    {
	/*
	 * For our European friends, A3 paper
	 * 
	 * For now, use (DEFAULT_BORDER * CMPERIN) as the best guess for a
	 * nonimageable area at the edges of the paper all around.
	 */

	dev_xmin = (DEFAULT_BORDER * CMPERIN) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_ymin = (DEFAULT_BORDER * CMPERIN) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_xmax = (42.0 - (DEFAULT_BORDER * CMPERIN)) * (1. / CMPERIN) * pixels_per_inch + .5;
	dev_ymax = (29.7 - (DEFAULT_BORDER * CMPERIN)) * (1. / CMPERIN) * pixels_per_inch + .5;
	ps_ypapersize = 29.7 * (1. / CMPERIN);
    }
    else
    {
	ps_xborder = -1.;
	ps_yborder = -1.;

	if (sscanf (ps_paper, "%fx%f%c,%f,%f",
		    &ps_ylength, &ps_xlength,
		    &units,
		    &ps_yborder, &ps_xborder) >= 2)
	{
	    if (units == 'i')
	    {
		/* Inches */
		unitscale = 1.;
	    }
	    else if (units == 'c')
	    {
		/* Centimeters */
		unitscale = 1. / CMPERIN;
	    }
	    else
	    {
		ERR (FATAL, name,
		     "Sorry, don't recognize units of type \"%c\"!", units);
	    }

	    /*
	     * Use a default border of DEFAULT_BORDER inches if they didn't
	     * specify it.
	     */
	    if (ps_xborder < 0.)
		ps_xborder = DEFAULT_BORDER / unitscale;

	    if (ps_yborder < 0.)
		ps_yborder = DEFAULT_BORDER / unitscale;


	    if (ps_xlength <= 0.)
	    {
		if (ps_ylength <= 0.)
		{
		    ERR (FATAL, name,
		    "Height and width of paper are both negative or zero!");
		}
		else
		    ps_xlength = 2. * ps_xborder + (ps_ylength - 2. * ps_yborder) / SCREEN_RATIO;
	    }
	    else if (ps_ylength <= 0.)
	    {
		ps_ylength = 2. * ps_yborder + (ps_xlength - 2. * ps_xborder) * SCREEN_RATIO;
	    }


	    ps_xlength *= unitscale;
	    ps_ylength *= unitscale;
	    ps_xborder *= unitscale;
	    ps_yborder *= unitscale;

	    ps_xmin = ps_xborder;
	    ps_ymin = ps_yborder;
	    ps_xmax = (ps_xlength - ps_xborder);
	    ps_ymax = (ps_ylength - ps_yborder);
	}
	else if (sscanf (ps_paper, "%f,%f,%fx%f,%f,%f%c",
			 &ps_ymin,
			 &ps_ymax,
			 &ps_ylength,
			 &ps_xmin,
			 &ps_xmax,
			 &ps_xlength,
			 &units) >= 6)
	{
	    if (units == 'i')
	    {
		/* Inches */
		unitscale = 1.;
	    }
	    else if (units == 'c')
	    {
		/* Centimeters */
		unitscale = 1. / CMPERIN;
	    }
	    else
	    {
		ERR (FATAL, name,
		     "Sorry, don't recognize units of type \"%c\"!", units);
	    }

	    ps_xmin *= unitscale;
	    ps_ymin *= unitscale;
	    ps_xmax *= unitscale;
	    ps_ymax *= unitscale;
	    ps_xlength *= unitscale;
	    ps_ylength *= unitscale;
	}
	else
	{
	    if (paperstring != NULL)
	    {
		ERR (COMMENT, name,
		     "Is your environmental variable \"$DEFAULT_PAPER_SIZE\" correct?");
	    }

	    ERR (FATAL, name,
		 "Don't recognize requested paper type \"%s\"!", ps_paper);
	}


	dev_xmin = ps_xmin * pixels_per_inch + .5;
	dev_ymin = ps_ymin * pixels_per_inch + .5;
	dev_xmax = ps_xmax * pixels_per_inch + .5;
	dev_ymax = ps_ymax * pixels_per_inch + .5;
	ps_ypapersize = ps_ylength;

	if (dev_xmin >= dev_xmax || dev_ymin >= dev_ymax)
	{
	    ERR (FATAL, name,
		 "Custom paper size has no plotting area.");
	}

	ps_set_papersize = YES;
    }


/*
 * device capabilities
 */
    getpar ("dumbfat", "1", &dumb_fat);
    getpar ("grey gray color", "1", &ps_color);
    getpar ("force", "1", &force_color);

    if (ps_color ||
	(strcmp (callname, "tcprpen") == 0) ||
	(strcmp (callname, "Tcprpen") == 0) ||
	(strcmp (wstype, "tcpr") == 0))
    {
	mono = NO;
	num_col = 16384;
    }
    else
    {
	mono = YES;
	num_col = 0;
    }
    getpar ("harddither", "1", &hard_dither);

    if (hard_dither)
    {
	smart_raster = YES;
	dev.raster = smart_psraster;
    }
    else
    {
	smart_raster = NO;
	dev.raster = dumb_psraster;
    }

    need_end_erase = YES;
    buffer_output = YES;
    smart_clip = NO;

    dither = 3;
    txfont = default_ps_font;
    txprec = DEFAULT_HARDCOPY_PREC;

    if (!GETPAR ("label", "s", label))
    {
	if ((user_name = getlogin ()) == NULL)
	    user_name = getpwuid (getuid ())->pw_name;
	dateget (date);
	sprintf (label, "%s, %s", user_name, date);
    }
    else
    {
	if ((strcmp (label, "no") == 0) || (strcmp (label, "n") == 0)
	    || (strcmp (label, " ") == 0) || (strcmp (label, "") == 0))
	    label[0] = '\0';
    }


/*
 *   Get the output file for the device
 */

    if (isatty (fileno (pltout)))
    {
/*
 * If it's a tty, then we certainly don't want to write to it!
 */
	creatafile = YES;
    }
    else
    {
/*
 *   If it's not a tty, they MAY be running pspen in a script
 *   and have logged out. In that case they want us to spool the plot,
 *   even though the output may not be a tty. (It may be /dev/null.)
 */
	if (fstat (fileno (pltout), &statbuf) != 0)
	{
	    ERR (WARN, name,
		 "Can't stat plot output! Trying to create a spool file instead.");

	    creatafile = YES;
	}
	else
	{
	    if (S_ISREG (statbuf.st_mode) || S_ISSOCK (statbuf.st_mode) ||
		S_ISFIFO (statbuf.st_mode))
	    {
		creatafile = NO;
	    }
	    else
	    {
		creatafile = YES;
	    }
	}
    }


    if (creatafile)
    {
char           *spooldirnm;

	file_created = YES;

	if ((spooldirnm = getenv ("VPLOTSPOOLDIR")) != NULL)
	{
	    sprintf (scratch_file, "%s%s", spooldirnm, "/PSPEN_XXXXXX");
	}
	else
	{
	    sprintf (scratch_file, "%s%s", PEN_SPOOL, "/PSPEN_XXXXXX");
	}

	mkstemp (scratch_file);
	if ((pltout = fopen (scratch_file, "w")) == NULL)
	{
	    if (spooldirnm != NULL)
	    {
		ERR (COMMENT, name,
		     "Check permissions on directory \"$VPLOTSPOOLDIR\" = \"%s\".",
		     spooldirnm);
	    }
	    else
	    {
		ERR (COMMENT, name,
		   "Perhaps you need to do \"setenv VPLOTSPOOLDIR /tmp\"?");
	    }

	    ERR (FATAL, name, "Could not open scratch file \"%s\"!",
		 scratch_file);
	}
    }

    if (GETPAR ("strip tex", "1", &yesget))
    {
	if (yesget)
	    tex = YES;
    }

    if (tex == YES)
    {
	/*
	 * allow as big a plot as we can and make sure our plot is in
	 * absolute coordinates. This prevents any clipping of the image. If
	 * tex=y is specified this plot will be scaled by TeX to fit whatever
	 * size we want, so don't clip it to page boundaries.
	 */
	dev_xmax = VP_MAX * RPERIN;
	dev_ymax = VP_MAX * RPERIN * SCREEN_RATIO;
	dev_xmin = 0;
	dev_ymin = 0;
	size = ABSOLUTE;
	label[0] = '\0';
    }

    if (GETPAR ("hold", "1", &yesget))
    {
	if (yesget)
	    hold = YES;
	else
	    hold = NO;
    }

    if (!GETPAR ("ncopies copies", "i", &ncopies_document))
    {
	ncopies_document = 1;
    }

/*
 * Initialize the PostScript file
 */
    if (tex == NO)
    {
	fprintf (pltout, "%%!\n");

	fprintf (pltout, "initgraphics 1 setlinecap 1 setlinejoin\n");
	fprintf (pltout, "%d rotate", 90);
	fprintf (pltout, " 0 %.2f translate %.2f %.2f scale gsave\n",
		 -ps_ypapersize * PSPERIN, psscale, psscale);
    }
    else
    {
	/*
	 * changed to work with psfig macros in TeX. Take the rotation and
	 * translation out. Now all we need is a shell that puts the
	 * BoundingBox stuff on the front.
	 */

	fprintf (pltout, " \n");
	fprintf (pltout, "%% Start of pspen output\n");
	fprintf (pltout, " \n");
	fprintf (pltout, " %.2f %.2f scale gsave\n", psscale, psscale);
	/*
	 * Make sure postscript uses the correct line join type to avoid
	 * "fanged S syndrome", though!
	 */
	fprintf (pltout, " 1 setlinecap 1 setlinejoin\n");
    }

    /* in case colorimage command is not known by this device */
    fprintf (pltout, " \n");
    fprintf (pltout, "systemdict /colorimage known not { \n");
    fprintf (pltout, "  /colortograyscale { \n");
    fprintf (pltout, "    dup /rgbdata exch store \n");
    fprintf (pltout, "    length 3 idiv \n");
    fprintf (pltout, "    /npixls exch store \n");
    fprintf (pltout, "    /indx 0 store \n");
    fprintf (pltout, "    /pixls npixls string store \n");
    fprintf (pltout, "    0 1 npixls -1 add { \n");
    fprintf (pltout, "      pixls exch \n");
    fprintf (pltout, "      rgbdata indx get .3 mul \n");
    fprintf (pltout, "      rgbdata indx 1 add get .59 mul add \n");
    fprintf (pltout, "      rgbdata indx 2 add get .11 mul add \n");
    fprintf (pltout, "      cvi \n");
    fprintf (pltout, "  put \n");
    fprintf (pltout, "      /indx indx 3 add store \n");
    fprintf (pltout, "    } for \n");
    fprintf (pltout, "    pixls \n");
    fprintf (pltout, "  } bind def \n");
    fprintf (pltout, "  /mergeprocs { \n");
    fprintf (pltout, "    dup length \n");
    fprintf (pltout, "    3 -1 roll \n");
    fprintf (pltout, "    dup \n");
    fprintf (pltout, "    length \n");
    fprintf (pltout, "    dup \n");
    fprintf (pltout, "    5 1 roll \n");
    fprintf (pltout, "    3 -1 roll \n");
    fprintf (pltout, "    add \n");
    fprintf (pltout, "    array cvx \n");
    fprintf (pltout, "    dup \n");
    fprintf (pltout, "    3 -1 roll \n");
    fprintf (pltout, "    0 exch \n");
    fprintf (pltout, "    putinterval \n");
    fprintf (pltout, "    dup \n");
    fprintf (pltout, "    4 2 roll \n");
    fprintf (pltout, "    putinterval \n");
    fprintf (pltout, "  } bind def \n");
    fprintf (pltout, "  /colorimage { \n");
    fprintf (pltout, "     pop \n");
    fprintf (pltout, "     pop \n");
    fprintf (pltout, "     {colortograyscale} \n");
    fprintf (pltout, "     mergeprocs \n");
    fprintf (pltout, "     image \n");
    fprintf (pltout, "  } bind def \n");
    fprintf (pltout, "} if \n");
    fprintf (pltout, " \n");

    /* Some abbreviations to make the files more compact */
    fprintf (pltout, "/m { moveto } bind def\n");
    fprintf (pltout, "/d { lineto } bind def\n");
    fprintf (pltout, "/r { rlineto } bind def\n");
    fprintf (pltout, "/x { 0 0 rlineto } bind def\n");
    fprintf (pltout, "/s { stroke } bind def\n");
    fprintf (pltout, "/gr { grestore } bind def\n");
    fprintf (pltout, "/np { newpath } bind def\n");
    fprintf (pltout, "/cf { closepath eofill } bind def\n");
    fprintf (pltout, " \n");

    if (!mono)
    {
	/*
	 * Begin by setting standard 7 colors. (The monochrome equivalent
	 * takes place in "psreset.c" for the monochrome case.)
	 */
	for (i = 0; i < 8; i++)
	{
	    psattributes (SET_COLOR_TABLE, i, red[i], green[i], blue[i]);
	}
	psattributes (SET_COLOR, 7, 0, 0, 0);
    }

    fflush (pltout);

    epause = 0;
    endpause = NO;
/*
 * This turns out not to be a good default for this hardcopy device,
 * since people think of it more as a screen-type device.
 *
 *  size = ABSOLUTE;
 */
}


/*
 * get the date
 */
#include <sys/time.h>
dateget (date)
    char           *date;
{
time_t          time ();
time_t          clock;
struct tm      *localtime ();
char           *asctime ();

    clock = time (0);
    sprintf (date, "%.16s", asctime (localtime (&clock)));
}
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
 *  source file:   ./filters/pslib/psplot.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), October 18 1988
 *	Stole imagplot.c to make psplot.c.
 * Joe Dellinger (SEP), October 19 1988
 *	Use "move/draw" method of sending paths.
 * Steve Cole (SEP), July 27 1989
 *      Added check on PATHLENGTH to avoid too many points/path.
 * Steve Cole (SEP), August 26 1990
 *      Removed check for staying at same point in addpath. Correct
 *	vplot behavior is to plot a point in this case, not ignore it.
 * Joe Dellinger (SOEST) June 23 1992
 *	The xold, yold variables in addpath were set but never used.
 * Joe Dellinger (Amoco) April 25 1996
 *	Use "s" instead of "stroke" to save a few bytes.
 *	Use "r" (short for "rlineto") when it uses less space than
 *	"d" (short for "lineto").
 * Joe Dellinger (Amoco) April 26 1996
 *	Jon has plots that consist of zillions of dots, so try to
 *	optimize that case a little bit better. ("x" is shorthand for "0 0 r")
 */

#include <stdio.h>
#include <string.h>
#include "../include/err.h"
#include "../include/enum.h"
#include "../include/extern.h"
#include "ps.h"

int             lost = 1;
static int      where = -1;
static int      ps_oldx = 0, ps_oldy = 0;

psplot (x, y, draw)
    int             x, y, draw;
{
    if (draw == 0)
    {
	startpath ();
	lost = 0;
    }
    addpath (x, y);
}

startpath ()
{
    if (where > 0)
    {
	endpath ();
    }
    where = 0;
}

addpath (x, y)
    int             x, y;
{
    /*
     * Just in case, allow lots of extra room
     */
char            stringd[80];
char            stringr[80];

    /* If where is -1 here it is a BUG! */
    if (where)
    {
	/*
	 * Important special case: a dot.
	 */
	if (x - ps_oldx == 0 && y - ps_oldy == 0)
	{
	    fprintf (pltout, "x\n");
	}
	else
	{
	    /*
	     * A bit brute force, but a sure way of using whichever is
	     * shorter: absolute or relative draws.
	     */
	    sprintf (stringd, "%d %d d\n", x, y);
	    sprintf (stringr, "%d %d r\n", x - ps_oldx, y - ps_oldy);

	    if (strlen (stringd) <= strlen (stringr))
		fprintf (pltout, "%s", stringd);
	    else
		fprintf (pltout, "%s", stringr);

	    fflush(pltout);
	}
    }
    else
	fprintf (pltout, "%d %d m\n", x, y);

    ps_oldx = x;
    ps_oldy = y;

    where++;
    if (where == PATHLENGTH - 1)
    {
	/* path reached maximum postscript pathlength; restart path */
	endpath ();
    }
}

endpath ()
{
    if (where > 0)
    {
	fprintf (pltout, "s\n");
	where = -1;
    }
    lost = 1;
}
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
 *  source file:   ./filters/pslib/psraster.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), September 11 1987
 *	Wrote dumb version of psraster.c.
 * Steve Cole (SEP), February 10 1988
 *	Rewrote using readhexstring instead of creating procedure for image.
 *      This avoids possible size problems.
 * Joe Dellinger, Feb 16 1988
 *	Make order and type of arguments same for both kinds of raster, to
 *	avoid problems with machines that have to have consistency of
 *	arguments.
 * Steve Cole (SEP), March 28 1988
 *      Cleaned up bit-handling code. Fixed bugs in handling of orient
 *      parameter and documented that code.
 * Steve Cole (SEP), March 31 1988
 *	Corrected polarity for an input 8 bits/pixel raster.
 * Joe Dellinger (SEP), October 19 1988
 *	Added "harddither" option.
 * Stew Levin (Mobil), May 8, 1996
 *	Changed mask0 definition to make source portable to machines
 *	that don't handle bit strings properly.
 */

#define mask0	((unsigned char) (((unsigned int) 1)  << 7))
#include <stdio.h>
#include "ps.h"
#include "../include/err.h"

extern FILE    *pltout;

dumb_psraster (count, out_of, xpos, ypos, length, orient, raster, dummy1, dummy2,byte2)
    int             count, out_of, xpos, ypos, length, orient, dummy1, dummy2,byte2;
    unsigned short  *raster;
{
int             i;
unsigned char   mask, outchar;
int             bitcount;
extern int      mono;
int             xshift, yshift;
int             xscale, yscale;
int             rangle;


/*   if(byte2==1)*/
/*   ERR(FATAL,"Asked for two byte option with dumb_psraster (unimplemented) \n");*/


    endpath ();

    if (count == 0)
    {
/*
 * Comments to explain the code for different orientations:
 * The PostScript command image outputs a raster with its lower left corner
 * at the origin, assuming that the raster is filled from (given the
 * transformation matrix we use below) left to right and top to bottom.
 * Depending on orient, this is not the way that the input raster is
 * filled. For orient 1 it is filled from top to bottom and right to
 * left, for 2 it is from right to left and bottom to top, and for 3
 * it is from left to right and bottom to top. We must apply a transformation
 * that changes the input raster so it appears to be filled in the orient
 * 0 manner. This is just a simple rotation. For orient=(0,1,2,3) the
 * rotation angle (measured clockwise) is (0,270,180,90) degrees.
 * This rotation is done about the origin.
 * To understand the postioning of the origin, for a given orientation,
 * start with a picture of what the raster looks like after the rotation.
 * The frist point of the first raster line (xpos,ypos) is at the upper
 * left. Now "back out" the rotation. The starting point is now at the
 * (upper left,upper right,lower right,lower left) corner for orient
 * (0,1,2,3). But we don't want the rotation to be performed about this
 * first point - we want it to be performed about the first point of the
 * LAST line of raster (sketch this on paper, then it is easy to see why).
 * So we want the origin to be at the beginning of the last line of raster.
 * Starting from the beginning of the first line, the shift needed to get
 * there is (0,-out_of) for orient 0, (-out_of,0) for orient 1,
 * (0,out_of) for orient 2, and (out_of,0) for orient 3. After this
 * translation and rotation, the scaling necessary to get the right
 * size raster is pretty obvious.
 *
 * Don't try to understand these comments unless you need to modify the
 * code! If you do need to modify the code, I think that these comments 
 * will help, even if it takes a few minutes to digest them.
 * A final note: the way to do all this neatly is to incorporate the
 * shifting, scaling, and rotating described here directly into the
 * transform matrix for image. But as long as this code works maybe
 * it's not worth the effort.
 * -S. Cole
 */

	if (orient == 0)
	{
	    xshift = xpos;
	    yshift = ypos - out_of;
	    xscale = length;
	    yscale = out_of;
	    rangle = 0;
	}
	else
	if (orient == 1)
	{
	    xshift = xpos - out_of;
	    yshift = ypos;
	    xscale = out_of;
	    yscale = length;
	    rangle = 270;
	}
	else
	if (orient == 2)
	{
	    xshift = xpos;
	    yshift = ypos + out_of;
	    xscale = length;
	    yscale = out_of;
	    rangle = 180;
	}
	else
	if (orient == 3)
	{
	    xshift = xpos + out_of;
	    yshift = ypos;
	    xscale = out_of;
	    yscale = length;
	    rangle = 90;
	}

	fprintf (pltout, "/picstr %d string def\n", (length + 7) / 8);

	fprintf (pltout, "gsave %d %d translate %d %d scale %d rotate\n", xshift, yshift, xscale, yscale, rangle);

	fprintf (pltout, "/raster {%d %d %d [ %d 0 0 %d 0 %d ] {currentfile picstr readhexstring pop} image} def\n", length, out_of, 1, length, -out_of, out_of);

	fprintf (pltout, "raster\n");
    }

    mask = mask0;
    outchar = 0;
    bitcount = 0;

    for (i = 0; i < length; i++)
    {
	if (raster[i] == 0)
	    outchar |= mask;
	mask >>= 1;
	bitcount++;
	if (bitcount == 8)
	{
	    fprintf (pltout, "%02.2x", outchar);
	    mask = mask0;
	    outchar = 0;
	    bitcount = 0;
	}
    }
    if (bitcount != 0)
	fprintf (pltout, "%02.2x", outchar);

    if (count != out_of - 1)
	return;

    fprintf (pltout, "\n");
    fprintf (pltout, "grestore\n");

}
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
 *  source file:   ./filters/pslib/psraster.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), September 11 1987
 *	Wrote smart version of psraster.c.
 * Steve Cole (SEP), February 10 1988
 *      Rewrote using readhexstring instead of creating procedure for image.
 *      This avoids possible size problems.
 *	Fixed orientation problems.
 * Steve Cole (SEP), March 28 1988
 *      Simplified the handling of the raster. In this routine we want
 *      always to output bytes and let PostScript do its own dithering.
 * Steve Cole (SEP), March 31 1988
 *      Corrected polarity of output raster.
 * Joe Dellinger (SEP), October 19 1988
 *	Added "harddither" option.
 * Steve Cole (SEP), November 30 1992
 *	Added color raster support.
 * Dave Nichols (SEP), April 7 1993
 *	Made grey rasters split their lines (like color ones do)
 * Joe Dellinger, David Lumley, 10-03-94
 * 	Check external variable ras_allgrey to decide if only gray raster
 *	is needed when color is requested by the pspen color=y option.
 */

#include <stdio.h>
#include "ps.h"
#include "../include/err.h"


extern FILE    *pltout;
extern int      ps_grey_ras[];
extern int	red[], green[], blue[];
extern int	mono;
extern int	ras_allgrey;

smart_psraster (xpix, ypix, xmin, ymin, xmax, ymax, raster_block, orient, dither_it,byte2)
    int             xpix, ypix, xmin, ymin, xmax, ymax, orient, dither_it,byte2;
    unsigned short  *raster_block;
{
int             i,j;
int             xshift, yshift;
int             xscale, yscale;
int             rangle;


/*   if(byte2==1)*/
/*   ERR(FATAL,"Asked for two byte option with smart_psraster (unimplemented) \n");*/


    endpath ();

    xscale = xmax - xmin;
    yscale = ymax - ymin;

    if (orient == 0)
    {
	xshift = xmin;
	yshift = ymin;
	rangle = 0;
    }
    else
    if (orient == 1)
    {
	xshift = xmin;
	yshift = ymax;
	rangle = 270;
    }
    else
    if (orient == 2)
    {
	xshift = xmax;
	yshift = ymax;
	rangle = 180;
    }
    else
    if (orient == 3)
    {
	xshift = xmax;
	yshift = ymin;
	rangle = 90;
    }

    fprintf (pltout, "gsave /picstr %d string def\n", xpix);

    fprintf (pltout, "%d %d translate %d %d scale %d rotate\n", xshift, yshift, xscale, yscale, rangle);

		fflush(pltout);
    if ( mono || ras_allgrey ) {
    fprintf (pltout, "/raster {%d %d 8 [ %d 0 0 %d 0 %d ] {currentfile picstr readhexstring pop} image} def\n", xpix, ypix, xpix, -ypix, ypix);

    fprintf (pltout, "raster\n");

    if (dither_it)
    {
	for (j = 0; j < xpix*ypix; j+=80)
	{
	   	for (i=j; (i<j+80 && i<xpix*ypix); i++)
	   	{
	    fprintf (pltout, "%02.2x", 255 - (int) raster_block[i]);
	   	}
    	   	fprintf (pltout, "\n");
	}
    }
    else
    {
	for (j = 0; j < xpix*ypix; j+=80)
	{
	   	for (i=j; (i<j+80 && i<xpix*ypix); i++)
	   	{
	    fprintf (pltout,"%02.2x", 255 - ps_grey_ras[(int) raster_block[i]]);
	   	}
    	   	fprintf (pltout, "\n");
	}
    }
    } else {
    fprintf (pltout, "/colraster {%d %d 8 [ %d 0 0 %d 0 %d ] {currentfile picstr readhexstring pop } false 3 colorimage} def\n", xpix, ypix, xpix, -ypix, ypix);

    fprintf (pltout, "colraster\n");

	for (j = 0; j < xpix*ypix; j+=80)
	{
		for (i=j; (i<j+80 && i<xpix*ypix); i++)
		{
	    	fprintf (pltout, "%02.2x%02.2x%02.2x", red[(int) raster_block[i]],green[(int) raster_block[i]],blue[(int) raster_block[i]]);
		}
    		fprintf (pltout, "\n");
	}
    }

    fprintf (pltout, "grestore\n");
}
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
 *  source file:   ./filters/pslib/psreset.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SEP), October 19 1988
 *	Added "grey" option.
 * Joe Dellinger (BP Amoco) October 5 1999
 *	Grey-level stuff should only be done for mono=y case.
 */

#include <stdio.h>
#include "../include/attrcom.h"
#include "ps.h"

extern FILE    *pltout;
extern int      mono;

psreset ()
{
int             ii;

    if (mono)
    {
	for (ii = 0; ii < 8; ii++)
	{
	    ps_grey_map (ii);
	}

	psattributes (SET_COLOR, 7, 0, 0, 0);
    }
}
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
 *  source file:   ./filters/pslib/pstext.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), September 16 1987
 *      Font 0 is now used if an undefined font is requested.
 * Steve Cole (SEP), March 23 1988
 *      For vplot fonts, added return after gentext call.
 * Steve Cole (SEP), April 2 1988
 *      Changed scaling to agree with other routines.
 * Steve Cole (SEP), June 20 1988
 *      Removed unused array "instruction".
 * Joe Dellinger (SOEST), March 1 1993
 *	Pspen hardware fonts were coming out too small due to a bug
 *	in this routine that used to not matter. This whole routine is
 *	full of bugs, though, and should be re-written from scratch!
 *	For example it should NOT attempt to use yscale and xscale itself.
 *	That is certainly WRONG! I took that part out; better to have
 *	text always come out square than have it scaled wrongly.
 *	(Talking to postscript in terms of "orient" and "size" is pointless
 *	and awkward in any case.) The "DIFFERENT_COORDINATES" stuff was not
 *	working and not worth saving, so I made it always "true".
 */

#include	<math.h>
#include	<stdio.h>
#include	"../include/vplot.h"
#include	"../include/extern.h"
#include	"../include/enum.h"
#include	"../include/params.h"
#include	"../include/round.h"
#include	"./ps.h"
#include	"./psfonts.h"

extern int      default_ps_font;

static double   path_orient_dx, path_orient_dy;
static double   up_orient_dx, up_orient_dy;

pstext (string, pathx, pathy, upx, upy)
    char           *string;
    float           pathx, pathy, upx, upy;
{
double          fpathx, fpathy, fupx, fupy;
double          up, path;
int             size, orient;
double          yfact, xfact;
static char     last_size = 0, last_font;

    endpath ();

    if (*string == '\0')
	return;

    if (txfont < 100)
    {
	gentext (string, pathx, pathy, upx, upy);
	return;
    }
    else if (txfont > psmaxfont)
    {
	txfont = default_ps_font;
    }

/*
 * Set the inital parameters
 */
    fpathx = (double) pathx;
    fpathy = (double) pathy;
    fupx = (double) upx;
    fupy = (double) upy;

    path = sqrt ((double) (fpathx * fpathx + fpathy * fpathy));
    up = sqrt ((double) (fupx * fupx + fupy * fupy));

    path_orient_dx = fpathx / path;
    path_orient_dy = fpathy / path;
    up_orient_dx = fupx / up;
    up_orient_dy = fupy / up;

/*
 * Postscript manual says height of "700 units" for the default
 * "1000 unit high" font is a typical "cap height". Vplot scales
 * by cap height, so we have to rescale by (1000./700.).
 * But I find that 570 is a better number!
 */
    size = ROUND (path * 1000. / 570.);
    orient = ROUND (acos (path_orient_dx) * 180 / 3.141592654);
    if (pathy < 0)
	orient *= -1;

/*
 *   Set the font and size
 */
    if ((size != last_size) || (txfont != last_font))
    {
	fprintf (pltout, "/%s findfont %d scalefont setfont\n",
		 psfonts[txfont - 100], size);
	last_size = size;
	last_font = txfont;
    }

    /*
     * SAVE THE CURRENT GRAPHICS STATE, ROTATE AND/OR SCALE THE COORDINATE
     * SYSTEM BEFORE MOVING
     */
    fprintf (pltout, "gsave\n");

    fprintf (pltout, "%d  %d translate\n", xold, yold);
    fprintf (pltout, "%d rotate\n", orient);


/*
 *  SET THE PROPER ALIGNMENT FROM THE CALCULATED LENGTH OF THE
 *  TEXT STRING
 */
    fprintf (pltout, "(%s) stringwidth pop\n", string);

    switch (txalign.ver)
    {
    case TV_TOP:
	yfact = -0.81 * (float) size;	/* was -1.0 */
	break;
    case TV_CAP:
	yfact = -0.654 * (float) size;	/* was -0.8 */
	break;
    case TV_SYMBOL:		/* CHECK THIS!!! */
    case TV_HALF:
	yfact = -0.327 * (float) size;	/* was -0.45 */
	break;
    case TV_BOTTOM:
	yfact = .1666666667 * (float) size;	/* was - */
	break;
    case TV_BASE:
	yfact = 0.0;
	break;
    default:
	break;
    }

    switch (txalign.hor)
    {
    case TH_RIGHT:
	xfact = -1.;
	break;
	break;
    case TH_NORMAL:
    case TH_LEFT:
	xfact = 0.;
	break;
    case TH_CENTER:
    case TH_SYMBOL:
	xfact = -.5;
	break;
    default:
	break;
    }
    fprintf (pltout, "%.2f mul\n", xfact);
    fprintf (pltout, "%.2f m\n", yfact);

    fprintf (pltout, "(%s) show\n", string);

    fprintf (pltout, "grestore\n");
}
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
 *  source file:   ./filters/pslib/psvector.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Steve Cole (SEP), March 24 1988
 *      Corrected computation of linewidth.
 * Steve Cole (SEP), March 31 1988
 *	Added line dashing logic.
 * Joe Dellinger (SEP) October 17 1988
 *	Stole imagvector.c to create psvector.c that does things with
 *	paths.
 * Joe Dellinger (SEP), October 19 1988
 *	Fixed bugs, use "move/draw" method of sending paths.
 * Aritomo Shinozaki (Loomis Laboratory of Physics), February 25 1989
 *	The code which generates dash pattern lines was incorrect.
 * Joe Dellinger (SOEST) June 23 1992
 *	Added dumb_fat option
 * Joe Dellinger (SOEST) October 15 1992
 *	No need for extra space at start of "setlinewidth" command.
 *	If the dash pattern has been completely lost (because of a
 *	grestore) call ps_set_dash to recreate it.
 */

#include <stdio.h>
#include "../include/extern.h"
#include "../include/enum.h"
#include "ps.h"
#define NEW_FAT		(nfat != ps_last_fat)

extern FILE    *pltout;
int             ps_last_fat = -1;

psvector (x1, y1, x2, y2, nfat, dashon)
    int             x1, y1, x2, y2;
    int             nfat, dashon;
{
static int      xlst, ylst;

/* Negative fatness? Skip it. */
    if (nfat < 0)
	return;

/* Clipped away? Skip it. */
    if (clip (&x1, &y1, &x2, &y2))
	return;

/*
 * ALAS, some postscript devices cannot be trusted to do their own
 * fattening. So we provide an option to do it all in the vplot
 * driver.
 */
    if (dumb_fat)
    {
	/* Have to dash FIRST */
	if (dashon)
	{
	    dashvec (x1, y1, x2, y2, nfat, dashon);
	    return;
	}

	/* Then fatten the resulting line segments SECOND */
	if (nfat)		/* Recursively calls itself to make fat lines */
	{
	    fatvec (x1, y1, x2, y2, nfat, dashon);
	    return;
	}

	/*
	 * (Unfortunately this means that if we do fattening in software we
	 * also have to do dashing in software.)
	 */
    }

/*
 *--------------------------------------------------------------------------
 * If dumb_fat == YES, then if we get to this point we must be
 * dealing with thin undashed lines. Just let the default code handle it;
 * it'll keep track of everything properly and won't try to do
 * anything because it will never need to be done.
 *--------------------------------------------------------------------------
 */

/*
 * set the fatness for the path
 */
    if (NEW_FAT)
    {
	endpath ();
	fprintf (pltout, "%d setlinewidth\n", nfat);
	ps_last_fat = nfat;
    }

/*
 * Make sure line dashing style is set correctly. The dashing pattern
 * was set in psattr, but it may or may not be turned on because
 * vplot can draw solid lines at any time without warning. Also,
 * if a new clipping window has been defined the dash pattern will have
 * been completely forgotten. So we have five cases to deal with:
 *
 * dashon=0 and !ps_dash_pattern_set -> thin line is default; just draw
 * dashon>0 and  ps_dash_pattern_set -> ready to draw dashed; just draw
 * dashon=0 and  ps_dash_pattern_set -> save the dash pattern
 *				        set solid line
 *				        ps_dash_pattern_set = NO
 *				        draw
 * dashon>0 and !ps_dash_pattern_set ->
 *				if (ps_dash_pattern_exists) then
 *					reload saved dash pattern
 *				        ps_dash_pattern_set = YES
 *				        draw
 *				else
 *					recreate lost pattern (ps_set_dash)
 *				        (sets ps_dash_pattern_set = YES)
 *				        (sets ps_dash_pattern_exists = YES)
 *					draw
 *				endif
 */
    if (dashon == 0 && ps_dash_pattern_set)
    {
	endpath ();
	fprintf (pltout, "currentdash\n");
	fprintf (pltout, "/vplotoffset exch def\n");
	fprintf (pltout, "/vplotdash exch def\n");
	fprintf (pltout, "[] 0 setdash\n");
	ps_dash_pattern_set = NO;
    }
    else if (dashon != 0 && !ps_dash_pattern_set)
    {
	if (ps_dash_pattern_exists)
	{
	    endpath ();
	    fprintf (pltout, "vplotdash vplotoffset setdash\n");
	    ps_dash_pattern_set = YES;
	}
	else
	{
	    endpath ();
	    /* This will set ps_dash_pattern_set and ps_dash_pattern_exists */
	    ps_set_dash (dashon);
	}
    }

    if ((x1 != xlst) || (y1 != ylst) || lost == 1)
    {
	/* Make sure it is a move, not a draw */
	psplot (x1, y1, 0);
    }
    psplot (x2, y2, 1);
    xlst = x2;
    ylst = y2;
}
