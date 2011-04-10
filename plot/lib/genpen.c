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

#include <rsf.h>

#include "genpen.h"
#include "device.h"
#include "polygon.h"

void vp_gen_area (vp_device dev, int npts, struct vp_vertex *head)
{
    struct vp_vertex  *v;
    int firstpoint, i;


    firstpoint = 2;
    vp_xminclip (0, 0, &firstpoint, dev);
    firstpoint = 1;

    v = head;
    for (i = 0; i < npts; i++) {
	vp_xminclip (v->x, v->y, &firstpoint,dev);
	v++;
    }

    firstpoint = -1;	/* Means this was the last point! */
    vp_xminclip (0, 0, &firstpoint, dev);
}

#ifdef kjhgkjhg

/*
extern int      cur_color;
extern int      need_devcolor;
*/

genhatch (npts, numhatch, angle, hafat, hacol, haoff, hasiz, head)
    int             npts, numhatch;
    struct vertex  *head;
    float           angle;
    int            *hafat, *hacol, *haoff, *hasiz;
{
register int    x, y, i;
int             xstr, xend, ystr, yend, y1, y2, x1, x2, start;
int             xa[4], ya[4], xwmin_r, xwmax_r, ywmin_r, ywmax_r;
int             skip, which_time;
int             ncross;
int             vminx, vmaxx, vminy, vmaxy;
struct vertex  *xhead, *yhead, *v;
int            *crosses;
int             cur_color_save;

    cur_color_save = cur_color;

    v = head;
    for (i = 0; i < npts; i++)
    {
	poly_rot (angle, &v->x, &v->y);
	v++;
    }

    /*
     * allocate storage for scan line cross points 
     */
    crosses = (int *) malloc ((unsigned) npts * sizeof (int));

    /*
     * double link the vertices. (head) is set to the node with the maximum
     * x-value so that intersect() will not eliminate 'head' while casting
     * off vertices. 
     */
    vminx = head->x;
    vmaxx = head->x;
    vminy = head->y;
    vmaxy = head->y;
    xhead = head;
    yhead = head;

    v = head;
    for (i = 0; i < npts; i++)
    {
	if (v->x > vmaxx)
	{
	    vmaxx = v->x;
	    xhead = v;
	}
	if (v->y > vmaxy)
	{
	    vmaxy = v->y;
	    yhead = v;
	}
	if (v->x < vminx)
	    vminx = v->x;
	if (v->y < vminy)
	    vminy = v->y;
	v++;
    }

/*
 * Find a new window which contains the old, rotated window
 */
    xwmin_r = xa[0] = xa[1] = xwmin;
    xwmax_r = xa[2] = xa[3] = xwmax;
    ywmin_r = ya[0] = ya[3] = ywmin;
    ywmax_r = ya[1] = ya[2] = ywmax;
    for (i = 0; i < 4; i++)
    {
	poly_rot (angle, &xa[i], &ya[i]);
	if (xwmin_r > xa[i])
	    xwmin_r = xa[i];
	if (ywmin_r > ya[i])
	    ywmin_r = ya[i];
	if (xwmax_r < xa[i])
	    xwmax_r = xa[i];
	if (ywmax_r < ya[i])
	    ywmax_r = ya[i];
    }

    if (vmaxx > xwmax_r)
	vmaxx = xwmax_r;
    if (vminx < xwmin_r)
	vminx = xwmin_r;
    if (vmaxy > ywmax_r)
	vmaxy = ywmax_r;
    if (vminy < ywmin_r)
	vminy = ywmin_r;

    /* stretch polygon in y-direction */
    v = yhead;
    do
    {
	v->y = 2 * (v->y) + 1;
	v = v->next;
    } while (v != yhead);

    for (which_time = numhatch; which_time < 2 * numhatch; which_time++)
    {
	if (hasiz[which_time] > 0)
	{
	    if (cur_color != hacol[which_time] || need_devcolor)
	    {
		cur_color = hacol[which_time];
		dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		need_devcolor = NO;
	    }

	    skip = hasiz[which_time];
	    start = haoff[which_time] + skip * (int) (vminy / skip);
	    for (y = start; y <= vmaxy; y += skip)
	    {
		ncross = intersect (2 * y, crosses, yhead, 1);
		sort (crosses, ncross);
		for (i = 0; i < ncross; i += 2)
		{
		    xstr = crosses[i];
		    xend = crosses[i + 1];
		    y1 = y2 = y;
		    poly_rot (-angle, &xstr, &y1);
		    poly_rot (-angle, &xend, &y2);
		    dev.vector (xstr, y1, xend, y2, hafat[which_time], 0);
		}
	    }
	}
    }
    /* shrink in y */
    v = yhead;
    do
    {
	v->y = ((v->y - 1) / 2);
	v = v->next;
    } while (v != yhead);

    /*
     * expand in x 
     */
    v = xhead;
    do
    {
	v->x = 2 * v->x + 1;
	v = v->next;
    } while (v != xhead);

    for (which_time = 0; which_time < numhatch; which_time++)
    {
	if (hasiz[which_time] > 1)
	{
	    if (cur_color != hacol[which_time] || need_devcolor)
	    {
		cur_color = hacol[which_time];
		dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		need_devcolor = NO;
	    }

	    skip = hasiz[which_time];
	    start = haoff[which_time] + skip * (int) (vminx / skip);
	    for (x = start; x <= vmaxx; x += skip)
	    {
		ncross = intersect (2 * x, crosses, xhead, 0);
		sort (crosses, ncross);
		for (i = 0; i < ncross; i += 2)
		{
		    ystr = crosses[i];
		    yend = crosses[i + 1];
		    x1 = x2 = x;
		    poly_rot (-angle, &x1, &ystr);
		    poly_rot (-angle, &x2, &yend);
		    dev.vector (x1, ystr, x2, yend, hafat[which_time], 0);
		}
	    }

	}
    }
    /*
     * shrink in x 
     */
    v = xhead;
    do
    {
	v->x = ((v->x - 1) / 2);
	v = v->next;
    } while (v != xhead);

    free ((char *) crosses);

    if (cur_color != cur_color_save)
    {
	cur_color = cur_color_save;
	need_devcolor = YES;
    }
}

poly_rot (angle, x, y)
    float           angle;
    int            *x, *y;
{
int             temp;
    temp = (*x * cos (angle) - *y * sin (angle)) + .5;
    *y = (*x * sin (angle) + *y * cos (angle)) + .5;
    *x = temp;
}

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

geninteract (what, controltty, string)
    int             what;
    FILE           *controltty;
    char           *string;
{
    switch (what)
    {
    case INT_PAUSE:
    case INT_F_PAUSE:
    case INT_GET_STRING:
	if (controltty == NULL)
	{
	    ERR (WARN, name, "Sorry, can't read string from terminal.");
	}
	else
	{
	    fgets (string, 79, controltty);
	}
	break;
    default:
	break;
    }
    return 0;
}

/*
 *
 *  source file:   ./filters/genlib/genmarker.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 *
 * Joe Dellinger, Dec 7 1987
 *	Fixed a bug that caused markers >127 and (>6 and <19)
 *	not to be drawn.
 * Joe Dellinger Feb 28 1988
 *	Even if size is zero still plot a point.
 */

#include "../include/params.h"
#include <stdio.h>
#include <math.h>
#include <vplot.h>
#include "../include/enum.h"
#include "../include/extern.h"
#include "../include/round.h"

genmarker (npts, type, size, pvec)
    int             npts, type, size;
    int            *pvec;
{
int             savetxfont, savetxprec, savetxovly, savefat;
struct txalign  savealign;
extern float    fatmult_orig;
char            txbuf[10];
char           *txbuffer = txbuf;

    savealign.hor = txalign.hor;
    savealign.ver = txalign.ver;

    savetxfont = txfont;
    savetxprec = txprec;
    savetxovly = txovly;
    savefat = fat;

/*
 * If it's shrunk away to nothing, plot it as a point.
 * Gentext eats zero-size text.
 */
    if (size == 0)
	type = 0;

    fat = (fatmult_orig * size * SYMBFATRATIO);
    fat += fatmult * (float) fatbase;
    txalign.hor = TH_SYMBOL;
    txalign.ver = TV_SYMBOL;
    txovly = OVLY_NORMAL;

    if ((type > 32) && (type < 127))	/* is it a printing character? */
    {
	*txbuffer = (char) (type);
	*(txbuffer + 1) = '\0';
	text_marker (txbuffer, size, npts, pvec);
    }
    else
    if (type >= 127)		/* special non-ASCII character */
    {
	sprintf (txbuffer, "\\v%d \0", type);
	text_marker (txbuffer, size, npts, pvec);
    }
    else			/* 0 through 5 are pre-defined; 6 through 20
				 * reserved */
    {
	switch (type)
	{
	case 0:
	case 1:
	    while (npts--)
	    {
		dev.point (*pvec, *(pvec + 1));
		pvec += 2;
	    }
	    break;
	case 2:		/* '+' */
	    txfont = MATH;
	    *txbuffer = (char) 57;
	    *(txbuffer + 1) = '\0';
	    text_marker (txbuffer, size, npts, pvec);
	    break;
	case 3:		/* '*' */
	    txfont = MATH;
	    *txbuffer = (char) 33;
	    *(txbuffer + 1) = '\0';
	    text_marker (txbuffer, size, npts, pvec);
	    break;
	case 4:		/* circle */
	    txfont = MISC;
	    *txbuffer = (char) 105;
	    *(txbuffer + 1) = '\0';
	    text_marker (txbuffer, size, npts, pvec);
	    break;
	case 5:		/* 'X' */
	    txfont = MATH;
	    *txbuffer = (char) 60;
	    *(txbuffer + 1) = '\0';
	    text_marker (txbuffer, size, npts, pvec);
	    break;
	case 20:		/* square */
	    txfont = MISC;
	    *txbuffer = (char) 72;
	    *(txbuffer + 1) = '\0';
	    text_marker (txbuffer, size, npts, pvec);
	    break;
	case 21:		/* triangle */
	    txfont = MISC;
	    *txbuffer = (char) 73;
	    *(txbuffer + 1) = '\0';
	    text_marker (txbuffer, size, npts, pvec);
	    break;
	case 22:		/* diamond */
	    txfont = MISC;
	    *txbuffer = (char) 74;
	    *(txbuffer + 1) = '\0';
	    text_marker (txbuffer, size, npts, pvec);
	    break;
	case 23:		/* star */
	    txfont = MISC;
	    *txbuffer = (char) 75;
	    *(txbuffer + 1) = '\0';
	    text_marker (txbuffer, size, npts, pvec);
	    break;
	default:
	    while (npts--)
	    {
		dev.point (*pvec, *(pvec + 1));
		pvec += 2;
	    }
	    break;
	}
    }

    txalign.hor = savealign.hor;
    txalign.ver = savealign.ver;
    txfont = savetxfont;
    txprec = savetxprec;
    txovly = savetxovly;
    fat = savefat;

    return (0);
}

text_marker (txbuffer, size, npts, pvec)
    char           *txbuffer;
    int             size, npts;
    int            *pvec;
{

    while (npts--)
    {
	xold = *pvec;
	yold = *(pvec + 1);
	if (txfont < NUMGENFONT)
	{
	    gentext (txbuffer,
	    /* Character path direction */
		     (float) size * aspect_ratio,
		     (float) 0,
	    /* Character up vector direction */
		     (float) 0,
		     (float) size);
	}
	else
	{
	    dev.text (txbuffer,
	    /* Character path direction */
		      (float) size * aspect_ratio,
		      (float) 0,
	    /* Character up vector direction */
		      (float) 0,
		      (float) size);
	}
	pvec += 2;
    }
}

/*
 *
 *  source file:   ./filters/genlib/genmessage.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

/*
 * Device independent subroutine to handle message operations
 */
#include	<stdio.h>
#include	"../include/mesgcom.h"
#include	"../include/enum.h"

genmessage (command, string)
    int             command;
    char            string[];
{

    switch (command)
    {
    case MESG_HOME:
    case MESG_READY:
	fflush (stderr);
	break;
    case MESG_TEXT:
	fprintf (stderr, "%s", string);
	break;
    case MESG_HIGHLIGHT_ON:
	/* Beep at them to get their attention */
	fprintf (stderr, "\07\07");
	break;
    case MESG_ON:
    case MESG_OFF:
    case MESG_ERASE:
    case MESG_HIGHLIGHT_OFF:
    case MESG_DONE:
    default:
	fflush (stderr);
	break;
    }
}

/*
 *
 *  source file:   ./filters/genlib/genpatarea.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger, Feb 16 1988
 *	Make number of arguments to dev.attributes and dev.raster consistent.
 */

#include <stdio.h>
#include "../include/pat.h"
#include "../include/vertex.h"
#include "../include/params.h"
#include "../include/extern.h"
#include "../include/enum.h"
#include "../include/attrcom.h"

/* Make the damn modulo function work for negative numbers */
#define MODULO(A,B)	(((A)>=0)?((A)%(B)):(((A)%(B)+(B))%(B)))

/*
 * This routine turns polygons filled with a pattern into calls to
 * dev.raster. (one line at a time).
 */
extern int      need_devcolor, cur_color;
extern char    *malloc ();

genpatarea (npts, head)
    int             npts;
    struct vertex  *head;
{
register int    y, i, ii;
register int    xstr, xend;
int             ncross;
int             vminx, vmaxx, vminy, vmaxy;
struct vertex  *yhead, *v;
int            *crosses;
unsigned char  *rasline;
static int      cur_color_save;
int             color;

/*
 * Save color so we can restore it the last time.
 */
    cur_color_save = cur_color;

    /*
     * allocate storage for scan line cross points 
     */
    crosses = (int *) malloc ((unsigned) npts * sizeof (int));

    /*
     * allocate storage for raster line 
     */
    rasline = (unsigned char *) malloc ((unsigned) (xwmax - xwmin + 1) * sizeof (unsigned char));

    /*
     * double link the vertices. (head) is set to the node with the maximum
     * x-value so that intersect() will not eliminate 'head' while casting
     * off vertices. 
     */
    vminx = head->x;
    vmaxx = head->x;
    vminy = head->y;
    vmaxy = head->y;
    yhead = head;

    v = head;
    for (i = 0; i < npts; i++)
    {
	if (v->x > vmaxx)
	{
	    vmaxx = v->x;
	}
	if (v->y > vmaxy)
	{
	    vmaxy = v->y;
	    yhead = v;
	}
	if (v->x < vminx)
	    vminx = v->x;
	if (v->y < vminy)
	    vminy = v->y;
	v++;
    }

    if (vmaxx > xwmax)
	vmaxx = xwmax;
    if (vminx < xwmin)
	vminx = xwmin;
    if (vmaxy > ywmax)
	vmaxy = ywmax;
    if (vminy < ywmin)
	vminy = ywmin;

    if ((pat[ipat] .ydim > 0) && (pat[ipat] .xdim > 0))
    {
	/* stretch polygon in y-direction */
	v = yhead;
	do
	{
	    v->y = 2 * (v->y) + 1;
	    v = v->next;
	} while (v != yhead);

	for (y = vminy; y <= vmaxy; y++)
	{
	    ncross = intersect (2 * y, crosses, yhead, 1);
	    sort (crosses, ncross);
	    for (i = 0; i < ncross; i += 2)
	    {
		xstr = crosses[i];
		xend = crosses[i + 1];

		if (xstr < xwmin && xend < xwmin)
		    continue;
		if (xstr > xwmax && xend > xwmax)
		    continue;

		if (xstr < xwmin)
		    xstr = xwmin;
		if (xend > xwmax)
		    xend = xwmax;

		if (pat[ipat] .xdim == 1)
		{
/* Faster to fill it with one vector */
		    color =
		     pat[ipat] .patbits[
					((pat[ipat] .ydim - 1) - (MODULO (y, pat[ipat] .ydim))) * pat[ipat] .xdim
		     ];
		    if (cur_color != color || need_devcolor)
		    {
			cur_color = color;
			dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
			need_devcolor = NO;
		    }
		    dev.vector (xstr, y, xend, y, 0, 0);
		}
		else
		{
		    for (ii = xstr; ii <= xend; ii++)
		    {
			rasline[ii - xstr] =
			 pat[ipat] .patbits[
					    ((pat[ipat] .ydim - 1) - (MODULO (y, pat[ipat] .ydim))) * pat[ipat] .xdim
					    + MODULO (ii, pat[ipat] .xdim)
			 ];
		    }
		    if (xstr <= xend)
			dev.raster (0, 1, xstr, y, xend - xstr + 1, 0, rasline, 0, 0);
		}
	    }
	}
	/* shrink in y */
	v = yhead;
	do
	{
	    v->y = ((v->y - 1) / 2);
	    v = v->next;
	} while (v != yhead);
    }

    free ((char *) rasline);
    free ((char *) crosses);

    if (cur_color != cur_color_save)
    {
	cur_color = cur_color_save;
	need_devcolor = YES;
    }
}

/*
 *
 *  source file:   ./filters/genlib/genpoint.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include "../include/extern.h"

genpoint (x1, y1)
    int             x1, y1;
{
    dev.vector (x1, y1, x1, y1, 0, 0);
}

/*
 *
 *  source file:   ./filters/genlib/genraster.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger Feb 16 1988
 *	Make number of arguments consistent between smart and dumb forms
 *	of dev.raster.
 */

#include <stdio.h>
#include "../include/extern.h"
#include "../include/enum.h"
#include "../include/attrcom.h"

extern int      overlay, cur_color, need_devcolor;

genraster (count, out_of, xpos, ypos, length, orient, raster, dummy1, dummy2)
    int             count, out_of, xpos, ypos, length, orient, dummy1, dummy2;
    unsigned char  *raster;
{
int             ii, sign, xy, xrpos, yrpos;
int             color, start;
static int      cur_color_save;

    switch (orient)
    {
    case 0:
	xrpos = xpos;
	yrpos = ypos;
	sign = 1;
	xy = 0;
	break;
    case 1:
	xrpos = ypos;
	yrpos = xpos;
	sign = -1;
	xy = 1;
	break;
    case 2:
	xrpos = xpos;
	yrpos = ypos;
	sign = -1;
	xy = 0;
	break;
    case 3:
	xrpos = ypos;
	yrpos = xpos;
	sign = 1;
	xy = 1;
	break;
    }

    start = xrpos;
    color = raster[0];

    if (count == 0)
    {
	/*
	 * First time remember the color so we can restore it the last time. 
	 */
	cur_color_save = cur_color;
    }

    for (ii = 0; ii <= length; ii++)
    {
	if (ii == length || raster[ii] != color)
	{
	    if (!((overlay == 1) && (color == 0)))
	    {
		if (cur_color != color || need_devcolor)
		{
		    cur_color = color;
		    dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		    need_devcolor = NO;
		}
		if (xy)
		    dev.vector (yrpos, start, yrpos, xrpos + sign * (ii - 1), 0, 0);
		else
		    dev.vector (start, yrpos, xrpos + sign * (ii - 1), yrpos, 0, 0);
	    }
	    if (ii != length)
	    {
		color = raster[ii];
		start = xrpos + sign * ii;
	    }
	}
    }

    if (count == out_of - 1 && cur_color != cur_color_save)
    {
	cur_color = cur_color_save;
	need_devcolor = YES;
    }
}

/*
 *
 *  source file:   ./filters/genlib/genraster1.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger Feb 16 1988
 *	Make number of arguments consistent between smart and dumb forms
 *	of dev.raster.
 */

#include <stdio.h>
#include "../include/extern.h"
#include "../include/enum.h"
#include "../include/attrcom.h"
#include "../include/params.h"

#define	MAXMEM 200000
#define OFF 0
#define ON 1

extern int      overlay, cur_color, need_devcolor, num_col_8;
extern char    *malloc ();
/*
 * Less than or equal to this, better to use point mode 
 */
int             break_point = 4;

/*
 * A more efficient genraster version for devices that can draw a
 * string of points quickly. Break the plot up by color, and by
 * points and vectors.
 */

genraster1 (count, out_of, xpos, ypos, length, orient, raster, dummy1, dummy2)
    int             count, out_of, xpos, ypos, length, orient, dummy1, dummy2;
    unsigned char  *raster;
{
int             ii, jj, yy, kk, ll;
int             xstart, state;
static int      cur_color_save;
static int      num_lines;
static unsigned char *array;
static int      ylength, ystart, jstart;
static int      color_used[MAX_COL + 1];
int             xsign, ysign, xy, xrpos, yrpos;

    switch (orient)
    {
    case 0:
	xrpos = xpos;
	yrpos = ypos;
	xsign = 1;
	ysign = 1;
	xy = 0;
	break;
    case 1:
	xrpos = ypos;
	yrpos = xpos;
	xsign = -1;
	ysign = 1;
	xy = 1;
	break;
    case 2:
	xrpos = xpos;
	yrpos = ypos;
	xsign = -1;
	ysign = -1;
	xy = 0;
	break;
    case 3:
	xrpos = ypos;
	yrpos = xpos;
	xsign = 1;
	ysign = -1;
	xy = 1;
	break;
    }

    if (count == 0)
    {
	/*
	 * First time remember the color so we can restore it the last time. 
	 */
	cur_color_save = cur_color;

	/*
	 * Also find out how many lines we can do at once. 
	 */
	num_lines = MAXMEM / length;
	if (num_lines < 1)
	    num_lines = 1;
	array = (unsigned char *) malloc ((unsigned) num_lines * length * sizeof (unsigned char));
	/*
	 * See whether we need to do color 0 or not 
	 */
	jstart = 0;
	if (overlay == 1)
	    jstart = 1;
    }

    if (ylength == 0)
    {
	/*
	 * Just starting a block. Remember where it is. 
	 */
	ystart = yrpos;
	for (ii = 0; ii <= num_col_8; ii++)
	{
	    color_used[ii] = 0;
	}
    }

/*
 * Save it.
 */
    for (ii = 0; ii < length; ii++)
    {
	array[length * ylength + ii] = raster[ii];
	color_used[raster[ii]] = 1;
    }
    ylength++;

    if (ylength >= num_lines || count == out_of - 1)
    {
/*
 * Plot it. Loop by color
 */

	need_devcolor = NO;

	for (kk = 0; kk < 2; kk++)
	{
/*
 * For maximum efficiency of output, better to have the kk loop
 * inside the jj loop. However, while watching on the screen better
 * to have "most important" stuff come out first.
 */
	    for (jj = jstart; jj < num_col_8; jj++)
	    {
		if (color_used[jj] == 0)
		    continue;

		cur_color = jj;
		dev.attributes (SET_COLOR, cur_color, 0, 0, 0);

		for (yy = 0; yy < ylength; yy++)
		{
		    state = OFF;
		    for (ii = 0; ii <= length; ii++)
		    {
			if (ii != length && array[length * yy + ii] == jj && state == OFF)
			{
			    xstart = xrpos + xsign * ii;
			    state = ON;
			    continue;
			}
			if ((ii == length || array[length * yy + ii] != jj) && state == ON)
			{
			    switch (kk)
			    {
			    case 0:
				if (xsign * (xrpos + xsign * ii - xstart) > break_point)
				{
				    if (xy)
					dev.vector (ystart - ysign * yy, xstart, ystart - ysign * yy, xrpos + xsign * (ii - 1), 0, 0);
				    else
					dev.vector (xstart, ystart - ysign * yy, xrpos + xsign * (ii - 1), ystart - ysign * yy, 0, 0);
				}
				break;
			    case 1:
				if (xsign * (xrpos + xsign * ii - xstart) <= break_point)
				{
				    for (ll = 0; ll < xsign * (xrpos + xsign * ii - xstart); ll++)
				    {
					if (xy)
					    dev.point (ystart - ysign * yy, xstart + xsign * ll);
					else
					    dev.point (xstart + xsign * ll, ystart - ysign * yy);
				    }
				}
				break;
			    }
			    state = OFF;
			}
		    }
		}

	    }
	}
	ylength = 0;
	if (count == out_of - 1)
	{
	    free ((char *) array);
	    if (cur_color != cur_color_save)
	    {
		cur_color = cur_color_save;
		need_devcolor = YES;
	    }
	}
    }
}

/*
 *
 *  source file:   ./filters/genlib/gentext.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */
/*
 * Joe Dellinger Oct 18 1987
 * 	Keep track of fatness as a float, not an int. Round when needed.
 * Joe Dellinger Jan 16 1988
 *	Allow user-defined fonts. As a check, require that they have a
 *	magic sequence on the front. If they ask for a font >= NUMGENFONT,
 *	modulo back into range.
 * Joe Dellinger Feb 16 1988
 *	Make number of arguments to dev.attributes consistent.
 * Joe Dellinger May 4 1988
 *	Explicitly check to make sure the ERROR font is compile-time
 *	loaded.
 * Joe Dellinger May 5 1990
 *	Shorts aren't worth the savings in file size for runtime-loaded
 *	fonts. Just make everything ints.
 * Joe Dellinger Oct 26 1991
 *	Fall back on looking for vplot binary font files in
 *	SYSTEM_FONT_DIRECTORY only if the environmental variable
 *	VPLOTFONTDIR is not set. Note both of these will normally
 *	be set to something ending in a "/"!
 *	Also, don't reject a file with a bad magic number out of
 *	hand, check to see if the problem is byte-swapping. If it
 *	is, just go ahead and silently fix it.
 */

/*
 * VPLOT soft text plotting
 *
 * Keywords: vplot text vector hershey font
 */


#include <sitedef.h>
#include	<stdio.h>
#include	<math.h>
#if defined(__stdc__) || defined(__STDC__)
#include <string.h>
#else
#include <strings.h>
#endif
#if defined(USG) || defined(SOLARIS)
#include       <fcntl.h>
#else /* USG */
#include	<sys/file.h>
#endif /* USG */
#include	<vplot.h>
#include	"../include/extern.h"
#include	"../include/err.h"
#include	"../include/enum.h"
#include	"../include/params.h"
#include	"../include/font_definitions.h"
#include	"../include/attrcom.h"
#include	"../include/round.h"

#define		NMARK   8	/* Maximum number of marks to use */
#define		MAXPOLY 100	/* Maximum number of points in polygon */

/*
 * The font to use for undefined fonts.
 * It should NOT be a runtime-loaded font!
 */
#define ERRFONT		0
/*
 * The glyph to use for undefined glyphs.
 * It must be a glyph in the font ERRFONT.
 * Needless to say, this glyph at least had better exist or you're
 * in real trouble.
 */
/* We use Glyph 30 in font 0, which is a '?' with a square around it. */
#define ERRGLYPH 	(30-font[0].dim[START])

/* Fraction of height of capital letter to use as vertical padding for box */
#define VSPACE_FRAC  .20
/* Fraction of inter-letter space to use for horizontal padding of box */
#define HSPACE_FRAC  0.5

#define EOC 0x8000	/* END OF CHARACTER BIT */
#define DRAWBIT 0x4000	/* DRAW BIT */
#define POLYBITS (EOC | DRAWBIT)	/* Polygon */
#define XBIT 0x0040
#define YBIT 0x2000
#define UNDEFINED	-1

#define SIGN_X (glyph_stroke & XBIT)
#define SIGN_Y (glyph_stroke & YBIT)
#define DRAW (glyph_stroke & DRAWBIT)
#define INPOLY ((glyph_stroke & POLYBITS) == POLYBITS)

#define COLOR_MAP(A) color_set[A][MAP];

/*
 * Correct for difference in size between what you want and what you've got
 */
#define SIZE_FACTOR(A) (((double)tsize/100.)*(A)/(double)(font[ttxfont].dim[CAP]-font[ttxfont].dim[BASE]))

/*
 * When you change sizes or fonts in midstream, what level do you line up on?
 */
#define ALIGN_HEIGHT	(font[ttxfont].dim[BASE])

#define CONTROL 001
#define BS	010
#define CR	015
#define NL	012

static double   path_orient_dx, path_orient_dy;
static double   up_orient_dx, up_orient_dy;
static double   xorigin_f, yorigin_f, xold_f, yold_f;
static int      ttxfont, cur_color_save, overlay_save;
extern int      cur_color, ipat, need_devcolor, overlay;
extern int      color_set[MAX_COL + 1][_NUM_PRIM];

#if defined(HAVE_STDLIB_H)
#include <stdlib.h>
#else
extern char    *malloc ();
extern char    *calloc ();
extern char    *getenv ();
#endif


static void
swab_font (stuff, bytecount)
    char           *stuff;
    int             bytecount;
{
int             icount;
char            temp;

    for (icount = 0;
	 icount < bytecount - sizeof (int) + 1;
	 icount += sizeof (int))
    {
	temp = stuff[icount + 0];
	stuff[icount + 0] = stuff[icount + 3];
	stuff[icount + 3] = temp;
	temp = stuff[icount + 1];
	stuff[icount + 1] = stuff[icount + 2];
	stuff[icount + 2] = temp;
    }
}

/*
 * interpret characters into vectors
 */
gentext (string, pathx, pathy, upx, upy)
    char           *string;
    float           pathx, pathy, upx, upy;
{
double          fpathx, fpathy, fupx, fupy;
double          up, path;
double          xp, yp;
float           tfat;
int             add;
int             ixp, iyp;
int            *istring;
unsigned int   *glyphptr, glyph_stroke;
double          xtxshift, ytxshift;
int             a, b, ii, jj, kk;
int             string_length;
double          last_widthl, last_widthr;
double          widthl_1st_char, char_width, total_width = 0.;
int             tsize, first, ttxfont_save, txfont_checked;
int             ghost = 0;
double          th_symb, th_symb_s, tv_symb, tv_symb_s;
double          char_width_s[NMARK];
double          total_width_s[NMARK];
int             mark_flag[NMARK];
double          last_widthl_s[NMARK], last_widthr_s[NMARK];
double          xold_f_s[NMARK];
double          yold_f_s[NMARK];
double          maxtop, minbot;
double          vspace, hspace, vline;
int             flag;
char           *charp;
int             polycount, xxx[MAXPOLY], yyy[MAXPOLY];
int            *ligp;
static int      one_error = YES;
int             linecount;

    if (*string == '\0')
	return;

/*
 * Set the initial parameters
 */

/*
 * Convert the input float path vectors to doubles.
 */
    fpathx = (double) pathx;
    fpathy = (double) pathy;
    fupx = (double) upx;
    fupy = (double) upy;

    path = sqrt ((double) (fpathx * fpathx + fpathy * fpathy));
    up = sqrt ((double) (fupx * fupx + fupy * fupy));

    if (path == 0. || up == 0.)
    {
/* Text got squashed away to nothing */
	return;
    }

    path_orient_dx = fpathx / path;
    path_orient_dy = fpathy / path;
    up_orient_dx = fupx / up;
    up_orient_dy = fupy / up;

/*
 * We didn't bomb out right away, so save things we may change so
 * they can be restored at exit
 */
    cur_color_save = cur_color;
    overlay_save = overlay;

    overlay = NO;

    for (ii = 0; ii < NMARK; ii++)
	mark_flag[ii] = 0;

/* Text starts out in the default font at 100% of the requested size */
    tsize = 100;

    if (font[ERRFONT].load == NO)
    {
	ERR (FATAL, name,
	  "(gentext) Font %d, the ERRFONT, must be loaded at compile time!",
	     ERRFONT);
    }

    txfont_checked = txfont;

    if (txfont_checked >= 0)
    {
	ttxfont = txfont_checked % NUMGENFONT;
	if (font[ttxfont].load == NO)
	{
	    load_font (txfont_checked);
	}
	txfont_checked = ttxfont;
    }
    else
    {
	txfont_checked = ERRFONT;
    }

    ttxfont = txfont_checked;
    total_width = 0.;
    char_width = font[ttxfont].dim[SPACE] * SIZE_FACTOR (path);
    th_symb = char_width / 2.;
    last_widthl = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
    last_widthr = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
    widthl_1st_char = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
    tv_symb = (font[ttxfont].dim[HALF] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
    maxtop = (font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
    minbot = (font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT) * SIZE_FACTOR (up);

/* These used in making the bounding box */
    vline = (font[ttxfont].dim[TOP] - font[ttxfont].dim[BOTTOM] + font[ttxfont].dim[LINE]) * SIZE_FACTOR (up);
    vspace = VSPACE_FRAC * (font[ttxfont].dim[CAP] - font[ttxfont].dim[BASE]) * SIZE_FACTOR (up);
    hspace = (HSPACE_FRAC * font[ttxfont].dim[SPACE] + font[ttxfont].dim[LETTER]) * SIZE_FACTOR (path);

/*
 * Parse ligatures, control sequences, etc.
 * each object gets 2 ints;
 * The first tells what sort of object it is.
 * Positive == printing character
 * Negative == command
 * The second is there for any parameters that are associated.
 * For normal characters, it gives the glyph number in this font.
 * For commands it gives any parameter the command may have.
 */
    string_length = strlen (string);
    istring = (int *) calloc ((unsigned) 2 * (string_length + 1), sizeof (int));

    for (ii = 0, charp = string; (charp - string) < string_length; ii++, charp++)
    {
	switch ((int) (*charp))
	{
/* Check for special ASCII characters first */
	case ' ':
	case NL:
	case CR:
	case BS:
	    istring[2 * ii] = -(int) (*charp);
	    continue;
	    break;
/* Check for \ commands */
	case '\\':
	    charp++;
	    switch (*charp)
	    {
/* \\ just falls through and makes a \ with no special properties */
	    case '\\':
		break;
/* \ commands with no arguments */
	    case '-':
	    case '>':
	    case '<':
	    case '^':
	    case '_':
	    case 'g':
	    case 'G':
	    case 'n':
	    case 'h':
		istring[2 * ii] = -(int) (*charp);
		continue;
		break;
/* \ commands with arguments */
	    case 's':
	    case 'f':
	    case 'F':
	    case 'k':
	    case 'r':
	    case 'm':
	    case 'M':
	    case 'v':
	    case 'c':
		istring[2 * ii] = -(int) (*charp);
		charp++;
/* default value of the argument is 0 if they just leave a space */
		istring[2 * ii + 1] = 0;
/* read the argument */
		sscanf (charp, "%d ", &istring[2 * ii + 1]);
/* skip past it and check for syntax */
		do
		{
		    if ((*charp >= '0' && *charp <= '9') ||
			*charp == '-' || *charp == '+')
		    {
			charp++;
		    }
		    else
			ERR (FATAL, name, "In text \\%c must be followed by an integer and then a space.", (char) (-istring[2 * ii]));
		} while (*charp != ' ');

		if (istring[2 * ii] == -(int) ('v'))
		{
/*
 * The \v command.
 * Make an ordinary character with the proper value.
 */
		    istring[2 * ii] = istring[2 * ii + 1];
		    istring[2 * ii + 1] = 0;
		}
		else
		if (istring[2 * ii] == -(int) ('F'))
		{
/* Font change command */
		    if (istring[2 * ii + 1] >= 0)
		    {
			ttxfont = istring[2 * ii + 1] % NUMGENFONT;
/* On this first pass through, load all the fonts we're going to need */
			if (font[ttxfont].load == NO)
			{
			    load_font (istring[2 * ii + 1]);
			}
			istring[2 * ii + 1] = ttxfont;
		    }
		    else
		    {
/* \F-1  means the default font again. */
			if (istring[2 * ii + 1] == -1)
			    istring[2 * ii + 1] = txfont_checked;
			else
			    istring[2 * ii + 1] = ERRFONT;
		    }
		}
		else
		if (istring[2 * ii] == -(int) ('c'))
		{
/* Color change command */
		    if (istring[2 * ii + 1] == -1)
		    {
/*
 * They want to return to the original text color.
 * This has already been checked to be within range and
 * properly mapped, so just use it!
 */
			istring[2 * ii + 1] = cur_color_save;
		    }
		    else
		    {
/*
 * Map from the color asked for to the colors that are available.
 * Normally only dovplot is allowed to do this, but this is an
 * unusual case where dovplot can't do it for us.
 */
			if (istring[2 * ii + 1] > MAX_COL || istring[2 * ii + 1] < 0)
			    ERR (FATAL, name, "(gentext) bad color number %d (max %d, min 0)",
				 istring[2 * ii + 1], MAX_COL);
			istring[2 * ii + 1] = COLOR_MAP (istring[2 * ii + 1]);
		    }
		}
		continue;
		break;
	    default:
		ERR (WARN, name, "(gentext) Unknown command \\%c.", *charp);
		charp--;
		break;
	    }
	default:
	    break;
	}
/* Normal character */
	istring[2 * ii] = (int) (*charp);
    }
    string_length = ii;

/* Ligatures */
    if (txprec > 1)
    {
	ttxfont = txfont_checked;
/*
 * Turning things into ligatures can only make the string shorter.
 * ii keeps track of where we are without ligatures,
 * kk keeps track of where we are with ligatures included.
 * The string is copied back into itself. Since ii >= kk, there is
 * no recursion problem.
 */
	for (ii = 0, kk = 0; ii < string_length; ii++, kk++)
	{
	    if (istring[2 * ii] < 0)
	    {
/*
 * The only special command we care about for constructing ligatures
 * is the font change command, since ligatures are font dependent.
 * The commands WILL break up a ligature, but other than that aren't
 * interpreted at all here.
 */
		if (-istring[2 * ii] == 'F')
		    ttxfont = istring[2 * ii + 1];
		istring[2 * kk] = istring[2 * ii];
		istring[2 * kk + 1] = istring[2 * ii + 1];
		continue;
	    }

/*
 * Take the first ligature that matches. This means that longer ligatures
 * MUST be listed first in the font data!
 */
/*
 * Loop over ligatures.
 * Each ligature has 1 number at the beginning giving the number of characters
 * in this ligature. The next number gives the glyph that is drawn for this
 * ligature. The next several numbers give the glyphs that are combined.
 * ligp points to the ligature we are currently searching for.
 * ligp += 2 + ligp[0] moves to the beginning of the next ligature:
 * 2 places for the 2 numbers at the beginning plus the ligp[0] characters
 * making up the ligature.
 */
	    for (ligp = font[ttxfont].lig; ligp[0] > 0; ligp += 2 + ligp[0])
	    {
/* Is there enough room before the end to possibly make this? */
		if (ii + ligp[0] - 1 < string_length)
		{
/* Loop over the characters in the ligature */
		    for (jj = 0; jj < ligp[0]; jj++)
		    {
/* Didn't match. Stop looking on this one. */
			if (ligp[jj + 2] != istring[2 * (ii + jj)])
			    goto failed;
		    }
/* Got to the end and so it worked. Put in the glyph for the ligature */
		    istring[2 * kk] = ligp[1];
/* skip past the ligp[0] characters in the original string that went into it */
		    ii += ligp[0] - 1;
		    goto success;
		}
	failed:
		continue;
	    }
/* No ligatures for this one. Copy it across unchanged */
	    istring[2 * kk] = istring[2 * ii];
/*
 * Don't need to look at any more ligatures for this character
 * (Ligatures don't nest)
 */
    success:
	    continue;
	}
/* Update the length of the string */
	string_length = kk;
    }


/************************************************************************/
/************************************************************************/

/*
 * This section conducts a "dry run" through the text string to determine
 * its length.
 *
 * Each character has a left half-width (widthl) and a right half-width (widthr).
 * These give the left and right half-widths of the character's bounding box
 * away from the character's origin.
 * The vertical dimensions of each character's bounding box are a function of
 * the font only; font[font_number].dim[TOP] and font[font_number].dim[BOTTOM]
 * give the coordinates in font units of the top and bottom of the character's box.
 *
 * Each character is also separated from its neighbors by an inter-letter space
 * (font[font_number].dim[LETTER]). This is effectively tacked onto the beginning
 * of each new character, with the exception of the first.
 *
 * When we actually output vectors we will start with the center of the
 * 1st character as our origin. (So we have to remember the left half-width
 * in order to be able to compensate for the fact that we measure the length
 * of the string starting from the left hand edge of the first character.)
 * Variables like "char_width" keep track of the latest character's width
 * in case we have to back up over it again.
 *
 * Multiplying by SIZE_FACTOR converts from the Font's units to Vplot's.
 * (Horizontal scales with "path", vertical scales with "up". These simply
 * give the length of the path and up vectors.)
 * All variables except for those in the font structure itself are in Vplot units.
 *
 * Left and right the total width of everything (characters and inter-spaces
 * between them) is summed into total_width. This is used to do the horizontal
 * text justification.
 *
 * Up and down the highest top and lowest bottom to date are saved in "maxtop" and
 * "minbot". These are used for vertical text justification. "ALIGN_HEIGHT"
 * gives the effective origin for vertical glyph positioning.
 *
 * The "symb" variables keep track of the symbol position of the latest character.
 */
    first = 1;
    flag = 1;
    linecount = 1;
    ttxfont = txfont_checked;
    for (ii = 0; ii < string_length; ii++)
    {
/* 
 * Figure the lenth of the message string.
 * Justification is based on the first line of text only
 * (That's what "flag" is for)
 */

/*
 * Check for special characters
 */
	if (istring[2 * ii] < 0)
	{
	    switch (-istring[2 * ii])
	    {
	    case 'n':
	    case NL:
		linecount++;
	    case CR:
		flag = 0;
		break;
	    case 'h':
	    case BS:
		total_width -= font[ttxfont].dim[LETTER] *
		 SIZE_FACTOR (path) + char_width;
		break;
	    case 'F':
/* Change the font */
		ttxfont = istring[2 * ii + 1];
		break;
	    case 's':
/* Change the size. This affects the SIZE_FACTOR. */
		tsize = istring[2 * ii + 1];
		break;
	    case 'k':
		if (flag)
		{
/*
 * Add in the blank space created by horizontal 'k'earning.
 * This is measured in percent of the width of a space in this font.
 *
 * Similar vertical movements are ignored for the purposes of justification.
 */
		    total_width += font[ttxfont].dim[SPACE]
		     * SIZE_FACTOR (path) * (istring[2 * ii + 1] / 100.);
		}
		break;
	    case 'm':
		if (istring[2 * ii + 1] < 0 || istring[2 * ii + 1] >= NMARK)
		    ERR (FATAL, name,
			 "(gentext) Too high a mark number %d", istring[2 * ii + 1]);
/* Save all relevant parameters as they are at this instant */
		if (flag)
		{
/* Vertical symbol alignment position */
		    tv_symb_s = tv_symb;
/* Horizontal symbol alignment position */
		    th_symb_s = th_symb;
/* Width of this character (in case the next thing is a backspace) */
		    char_width_s[istring[2 * ii + 1]] = char_width;
/* The width so far up to this point */
		    total_width_s[istring[2 * ii + 1]] = total_width;
		}
		mark_flag[istring[2 * ii + 1]] = 1;
		break;
	    case 'M':
		if (istring[2 * ii + 1] < 0 || istring[2 * ii + 1] >= NMARK)
		    ERR (FATAL, name,
			 "(gentext) Too high a mark number %d", istring[2 * ii + 1]);
/* Make sure it isn't junk */
		if (!mark_flag[istring[2 * ii + 1]])
		    ERR (FATAL, name,
			 "(gentext) Attempt to use undefined mark number %d",
			 istring[2 * ii + 1]);
/*
 * Restore the parameters previously saved. All events after that point
 * are now ignored for the purposes of justification.
 */
		if (flag)
		{
		    tv_symb = tv_symb_s;
		    th_symb = th_symb_s;
		    char_width = char_width_s[istring[2 * ii + 1]];
		    total_width = total_width_s[istring[2 * ii + 1]];
		}
		break;
	    case '-':		/* Nothing */
		break;
	    case '>':		/* Forward one inter-letter space */
		if (flag)
		    total_width += font[ttxfont].dim[LETTER]
		     * SIZE_FACTOR (path);
		break;
	    case '<':		/* Remove one inter-letter space */
		if (flag)
		    total_width -= font[ttxfont].dim[LETTER]
		     * SIZE_FACTOR (path);
		break;
	    case '^':		/* Up a half letter */
	    case '_':		/* Down a half letter */
	    case 'g':		/* Make text invisible */
	    case 'G':		/* Make text visible again */
		break;
	    case ' ':		/* Space */
		if (flag)
		{
		    char_width = font[ttxfont].dim[SPACE]
		     * SIZE_FACTOR (path);
		    th_symb = char_width / 2.;
		    tv_symb = (font[ttxfont].dim[HALF] - ALIGN_HEIGHT)
		     * SIZE_FACTOR (up);
		    total_width += char_width;
		    if (first)
		    {
/* If it is the first character, remember the left half width */
			widthl_1st_char = font[ttxfont].dim[SPACE] * .5
			 * SIZE_FACTOR (path);
/* No longer at the first printable character */
			first = 0;
		    }
		    else
		    {
/* else add inter-letter space between it and the previous character */
			total_width += font[ttxfont].dim[LETTER]
			 * SIZE_FACTOR (path);
		    }
		}
		break;
	    }
	    continue;
	}

	if (flag)
	{
/*
 * There are 2 ways a glyph can be undefined: it can be outside the range of
 * the font, OR it can have no data associated with it
 */
	    if (istring[2 * ii] >= font[ttxfont].dim[START] && istring[2 * ii] <= font[ttxfont].dim[END])
	    {
/* Find the glyph number, and save it so we don't have to recalculate it */
		istring[2 * ii + 1] = istring[2 * ii] - font[ttxfont].dim[START];
		if (font[ttxfont].saddr[istring[2 * ii + 1]] != UNDEFINED)
		{
/* OK glyph */
/* In case it's the last one, save its vertical symbol position */
		    tv_symb = (font[ttxfont].symbol[istring[2 * ii + 1]] - ALIGN_HEIGHT)
		     * SIZE_FACTOR (up);
/* And in case we back up later its width */
		    char_width = (font[ttxfont].swidthl[istring[2 * ii + 1]] +
				  font[ttxfont].swidthr[istring[2 * ii + 1]])
		     * SIZE_FACTOR (path);
/* and horizontal symbol position */
		    th_symb = font[ttxfont].swidthr[istring[2 * ii + 1]]
		     * SIZE_FACTOR (path);
/* See if it sets a new record high */
		    if ((font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up) > maxtop)
			maxtop = (font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
/* Or a record low */
		    if ((font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT)
			* SIZE_FACTOR (up) < minbot)
			minbot = (font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT)
			 * SIZE_FACTOR (up);
/* Add it into the total width */
		    total_width += char_width;
		    if (first)
		    {
/* If it's the first remember its left half width */
			widthl_1st_char = font[ttxfont].swidthl[istring[2 * ii + 1]]
			 * SIZE_FACTOR (path);
		    }
		    else
		    {
/* or if not first add in the space between it and the previous glyph */
			total_width += font[ttxfont].dim[LETTER]
			 * SIZE_FACTOR (path);
		    }
		}
		else
		{
/* Second way to be undefined. Turn it into a "special" character */
		    istring[2 * ii] = UNDEFINED;
		}
	    }
	    else
	    {
/* First way to be undefined. Turn it into a "special" character */
		istring[2 * ii] = UNDEFINED;
	    }

	    if (istring[2 * ii] == UNDEFINED)
	    {
/*
 * If it is undefined, use the special "ERROR" glyph and then
 * treat that just like we would treat a regular character.
 */
		ttxfont_save = ttxfont;
		ttxfont = ERRFONT;
		istring[2 * ii + 1] = ERRGLYPH;
		char_width = (font[ttxfont].swidthl[ERRGLYPH] +
			      font[ttxfont].swidthr[ERRGLYPH])
		 * SIZE_FACTOR (path);
		th_symb = font[ttxfont].swidthr[ERRGLYPH] * SIZE_FACTOR (path);
		tv_symb = (font[ttxfont].symbol[ERRGLYPH] - ALIGN_HEIGHT)
		 * SIZE_FACTOR (up);
		if ((font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up) > maxtop)
		    maxtop = (font[ttxfont].dim[TOP] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
		if ((font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT) * SIZE_FACTOR (up) < minbot)
		    minbot = (font[ttxfont].dim[BOTTOM] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
		total_width += char_width;
		if (first)
		{
		    widthl_1st_char = font[ttxfont].swidthl[ERRGLYPH] * SIZE_FACTOR (path);
		}
		else
		    total_width += font[ttxfont].dim[LETTER] * SIZE_FACTOR (path);
		ttxfont = ttxfont_save;
	    }

/* We printed something, so we aren't at the first character anymore */
	    first = 0;
	}
	else
	{
/*
 * If we're past the first line of text, do the few things that aren't related
 * to justification (looking for undefined glyphs, finding the glyph numbers)
 */
	    if (istring[2 * ii] >= font[ttxfont].dim[START] && istring[2 * ii] <= font[ttxfont].dim[END])
	    {
		istring[2 * ii + 1] = istring[2 * ii] - font[ttxfont].dim[START];
		if (font[ttxfont].saddr[istring[2 * ii + 1]] == UNDEFINED)
		{
		    istring[2 * ii] = UNDEFINED;
		}
	    }
	    else
	    {
		istring[2 * ii] = UNDEFINED;
	    }
	    if (istring[2 * ii] == UNDEFINED)
	    {
		istring[2 * ii + 1] = ERRGLYPH;
	    }
	}
    }

/*
 *  Set the proper alignment from the calculated length of the 
 *  text string. Remember that when we plot zero will be in the center
 *  of the first character and not the left hand edge of the first character,
 *  so we have to use widthl_1st_char to compensate for that.
 */
    switch (txalign.hor)
    {
    case TH_SYMBOL:
	xtxshift = total_width - widthl_1st_char - th_symb;
	break;
    case TH_CENTER:
	xtxshift = total_width / 2. - widthl_1st_char;
	break;
    case TH_RIGHT:
	xtxshift = total_width - widthl_1st_char;
	break;
    case TH_NORMAL:
    case TH_LEFT:
    default:
	xtxshift = -widthl_1st_char;
	break;
    }

    tsize = 100;
    ttxfont = txfont_checked;
    tfat = fat;

/*
 * CAP, HALF, and BASE are calculated based on font and size of default
 * TOP and BOTTOM are based on highest TOP and lowest BOTTOM of all
 * glyphs in string.
 */
    switch (txalign.ver)
    {
    case TV_SYMBOL:
	ytxshift = tv_symb;
	break;
    case TV_TOP:
	ytxshift = maxtop;
	break;
    case TV_CAP:
	ytxshift = (font[ttxfont].dim[CAP] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
	break;
    case TV_HALF:
	ytxshift = (font[ttxfont].dim[HALF] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
	break;
    case TV_BOTTOM:
	ytxshift = minbot;
	break;
    case TV_NORMAL:
    case TV_BASE:
    default:
	ytxshift = (font[ttxfont].dim[BASE] - ALIGN_HEIGHT) * SIZE_FACTOR (up);
	break;
    }


/************************************************************************/
/************************************************************************/


/*
 * This part of the code draws the characters.
 *
 * The complexity arises because when we do each character we have to
 * be at its CENTER in the left-right direction and at its ALIGN_HEIGHT
 * in the up-down direction.
 * So to move from one character to the next we have to move over by
 * the right half-width of the previous character, an inter-letter space,
 * and then the left half-width of the character we're on!
 * This helps to make the code confusing.
 *
 * We also have to always be ready to unexpectedly back up over the last
 * character, so we also have to keep around the left half-width of the
 * previous character as well.
 */
    xold_f = xold;
    yold_f = yold;

/*
 * The "mov" routine moves us around in the cock-eyed coordinate system
 * determined by the up and path vectors. There's no problem if these
 * vectors aren't orthogonal!
 */
    mov (-xtxshift, -ytxshift);
    xorigin_f = xold_f;
    yorigin_f = yold_f;

    if (txovly)
    {
	mov (-widthl_1st_char - hspace, minbot - vline * (linecount - 1) - vspace);
	xxx[0] = ROUND (xold_f);
	yyy[0] = ROUND (yold_f);
	mov (0., maxtop - minbot + vline * (linecount - 1) + 2. * vspace);
	xxx[1] = ROUND (xold_f);
	yyy[1] = ROUND (yold_f);
	mov (total_width + 2. * hspace, 0.);
	xxx[2] = ROUND (xold_f);
	yyy[2] = ROUND (yold_f);
	mov (0., -(maxtop - minbot + vline * (linecount - 1) + 2. * vspace));
	xxx[3] = ROUND (xold_f);
	yyy[3] = ROUND (yold_f);

	if (txovly == 2 || txovly == 3)
	{
	    if (cur_color != 0 || need_devcolor)
	    {
		cur_color = 0;
		dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		need_devcolor = NO;
	    }
	    drawpolygon (4, xxx, yyy);
	    if (cur_color != cur_color_save || need_devcolor)
	    {
		cur_color = cur_color_save;
		dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		need_devcolor = NO;
	    }
	}

	if (txovly == 1 || txovly == 3)
	{
	    dev.vector (xxx[0], yyy[0], xxx[1], yyy[1], ROUND (tfat), 0);
	    dev.vector (xxx[1], yyy[1], xxx[2], yyy[2], ROUND (tfat), 0);
	    dev.vector (xxx[2], yyy[2], xxx[3], yyy[3], ROUND (tfat), 0);
	    dev.vector (xxx[3], yyy[3], xxx[0], yyy[0], ROUND (tfat), 0);
	}

	xold_f = xorigin_f;
	yold_f = yorigin_f;
    }

    first = 1;

/*
 * This is where the actual drawing of the characters takes place.
 * Loop over all the characters in the string
 */
    for (ii = 0; ii < string_length; ii++)
    {
/*
 * Check for special characters first
 */
	if (istring[2 * ii] < 0)
	{
	    switch (-istring[2 * ii])
	    {			/* standard carriage controls */
	    case 'h':
	    case BS:
		mov (-font[ttxfont].dim[LETTER] * SIZE_FACTOR (path) -
		     (last_widthl + last_widthr), 0.);
		break;
	    case NL:
	    case 'n':
		xold_f = xorigin_f;
		yold_f = yorigin_f;
		mov (0., -(
			   (font[ttxfont].dim[TOP] - font[ttxfont].dim[BOTTOM] + font[ttxfont].dim[LINE])
			   * SIZE_FACTOR (up)));
		xorigin_f = xold_f;
		yorigin_f = yold_f;
		first = 1;
		break;
	    case CR:
		xold_f = xorigin_f;
		yold_f = yorigin_f;
		first = 1;
		break;
	    case 'F':
		ttxfont = istring[2 * ii + 1];
		break;
	    case 'c':
		if (cur_color != istring[2 * ii + 1] || need_devcolor)
		{
		    cur_color = istring[2 * ii + 1];
		    dev.attributes (SET_COLOR, cur_color, 0, 0, 0);
		    need_devcolor = NO;
		}
		break;
	    case 's':
		tsize = istring[2 * ii + 1];
		break;
	    case 'f':
		tfat += istring[2 * ii + 1] * fatmult;
		break;
	    case 'k':
/* Horizontal motion */
		mov (font[ttxfont].dim[SPACE] * SIZE_FACTOR (path) *
		     (istring[2 * ii + 1] / 100.), 0.);
		break;
	    case 'r':
/* Vertical motion */
		mov (0., (font[ttxfont].dim[CAP] - font[ttxfont].dim[BASE])
		     * SIZE_FACTOR (up) * (istring[2 * ii + 1] / 100.));
		break;
	    case 'm':
/*
 * Save the current position. No need to check mark number valid; this
 * was checked when we did the justification
 */
		last_widthl_s[istring[2 * ii + 1]] = last_widthl;
		last_widthr_s[istring[2 * ii + 1]] = last_widthr;
		xold_f_s[istring[2 * ii + 1]] = xold_f;
		yold_f_s[istring[2 * ii + 1]] = yold_f;
		break;
	    case 'M':
/*
 * Restore the current position
 */
		last_widthl = last_widthl_s[istring[2 * ii + 1]];
		last_widthr = last_widthr_s[istring[2 * ii + 1]];
		xold_f = xold_f_s[istring[2 * ii + 1]];
		yold_f = yold_f_s[istring[2 * ii + 1]];
		break;
	    case 'G':
		ghost = 0;
		break;
	    case 'g':
		ghost = 1;
		break;
	    case '^':
/* Up half a character */
		mov (0.,
		     (font[ttxfont].dim[CAP] - font[ttxfont].dim[BASE])
		     * (.5) * SIZE_FACTOR (up));
		break;
	    case '_':
/* Down half a character */
		mov (0., -(
			   (font[ttxfont].dim[CAP] - font[ttxfont].dim[BASE])
			   * (.5) * SIZE_FACTOR (up)));
		break;
	    case '-':
		break;
	    case '>':
/* Right an inter-letter space */
		mov (font[ttxfont].dim[LETTER] * SIZE_FACTOR (path), 0.);
		break;
	    case '<':
/* Left an inter-letter space */
		mov (-(font[ttxfont].dim[LETTER] * SIZE_FACTOR (path)), 0.);
		break;
	    case -(UNDEFINED):
		/* Don't overload them with error messages */
		if (one_error)
		{
		    ERR (WARN, name,
			 "(gentext) Attempt(s) to use undefined glyph(s) in font %d",
			 ttxfont);
		    one_error = NO;
		}
/* Switch to use the ERROR glyph, and the treat it as a regular glyph */
		ttxfont_save = ttxfont;
		ttxfont = ERRFONT;
		goto not_special;
		break;
	    case ' ':
	    default:
		if (!first)
		{
		    mov (
			 (font[ttxfont].dim[SPACE] * .5 +
			  font[ttxfont].dim[LETTER])
			 * SIZE_FACTOR (path)
			 + last_widthr
			 ,0.);
		}
		else
		{
		    first = 0;
		}
		last_widthl = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
		last_widthr = font[ttxfont].dim[SPACE] * .5 * SIZE_FACTOR (path);
		break;
	    }
	}
	else
	{
    not_special:
/*
 *  Printable character.
 *  Pull out the actual strokes that make up each glyph.
 *  First get the address of the character from the address array
 *  Get the address by adding the offset (add) to the base array address
 *  (font[ttxfont].svec).
 */
	    add = font[ttxfont].saddr[istring[2 * ii + 1]];
	    glyphptr = font[ttxfont].svec + add;
/*
 * Now that we have the address of the fonts,
 * we position the pen at the beginning of the next character
 * (Unless it's the first printable character in which case we are already
 *  there)
 */
	    if (!first)
	    {
		mov (
		     (font[ttxfont].swidthl[istring[2 * ii + 1]] +
		      font[ttxfont].dim[LETTER])
		     * SIZE_FACTOR (path)
		     + last_widthr
		     ,0.);
	    }
	    else
		first = 0;

/* Save the left and right half-widths of this glyph */
	    last_widthl = font[ttxfont].swidthl[istring[2 * ii + 1]]
	     * SIZE_FACTOR (path);
	    last_widthr = font[ttxfont].swidthr[istring[2 * ii + 1]]
	     * SIZE_FACTOR (path);

/*
 * Calculate where to position each character in high precision.
 */
	    xnew = ROUND (xold_f);
	    ynew = ROUND (yold_f);
/*
 * This loop contains the structure for the actual drawing of the characters
 * We go through this block until an "END OF CHARACTER" is read
 *
 *  Strokes are kept in a packed format to save
 *  space. Each stroke is packed into an unsigned
 *  short int with the following format: 
 *
 *         edsyyyyyysxxxxxx
 *         ||||____|||____|--> The x-coordinate value
 *         |||   |  `--------> Set if X < 0			
 *         |||   `-----------> The y-coordinate value
 *         ||`---------------> Set if Y < 0
 *         |`----------------> Draw bit, set if
 *         |                    command is draw.
 *         |                    Clear, if move.
 *         `-----------------> End of Character Bit
 *
 *  This is enough bits per coordinate to accomodate all the "Hershey" fonts.
 *
 *  Polygons are also encoded into this scheme. If the EOC and DRAW bits
 *  are simultaneously on, then we are inside a polygon. The last point of
 *  the polygon will only have the draw flag on, and at that point the entire
 *  polygon will be outputted.
 *
 *  Unfortunately it was later discovered that using shorts and ints mixed
 *  together was a bad idea, so everything is now ints even though the
 *  top part of some ints are unused!!! So we get the worst of both worlds.
 */

	    polycount = 0;

	    while ((glyph_stroke = *glyphptr++) != EOC)
	    {
		a = glyph_stroke & 077;
		if (SIGN_X)
		    a = -a;
		b = (glyph_stroke >> 7) & 077;
		if (SIGN_Y)
		    b = -b;
		b -= ALIGN_HEIGHT;

/*
 * Here is the correct place to insert code to rotate a glyph.
 * You want to do that before it gets distorted by the global coordinate
 * transformation. The "ALIGN_HEIGHT" defines where the vertical origin
 * is for a glyph in a font. Note that we are in the font's coordinate
 * system units at this point.
 */
		xp = xold_f;
		yp = yold_f;
/*
 * Cock-eyed coordinate system.
 * "up" is in the direction of the up vector,
 * "right" is in the direction of the path vector.
 * These can be screwy, and thus distort the glyph.
 *
 * "path" and "up" contain the magnitude, and the
 * "orient_dx"'s and "orient_dy"'s contain the direction cosines
 */
		xp += SIZE_FACTOR (path) * a * path_orient_dx +
		 SIZE_FACTOR (up) * b * up_orient_dx;
		yp += SIZE_FACTOR (path) * a * path_orient_dy +
		 SIZE_FACTOR (up) * b * up_orient_dy;
		ixp = ROUND (xp);
		iyp = ROUND (yp);

		if (polycount > 0 && !INPOLY)
		{
/*
 * If we just WERE in a polygon, but are not now, then we must have just
 * finished one. Plot it out.
 */
		    xxx[polycount] = ixp;
		    yyy[polycount] = iyp;
		    polycount++;
		    if (!ghost)
		    {
			drawpolygon (polycount, xxx, yyy);
		    }
/*
 * Done with this one, reset the vertex counter
 */
		    polycount = 0;
		}
		else
		if (INPOLY)
		{
/*
 * We're still saving up the polygon. Save this vertex.
 */
		    xxx[polycount] = ixp;
		    yyy[polycount] = iyp;
		    polycount++;
		    if (polycount > MAXPOLY - 1)
		    {
			ERR (FATAL, name,
			     "(gentext) Too many points in polygon.");
		    }
		}
		else
		{
/* No polygons, just output the vector */
		    if (DRAW && !ghost)
			dev.vector (xnew, ynew, ixp, iyp, ROUND (tfat), 0);
		}

		xnew = ixp;
		ynew = iyp;
	    }

/*
 * If we switched to the error font just to plot the error glyph,
 * then switch back to the correct font
 */
	    if (istring[2 * ii] == UNDEFINED)
	    {
		ttxfont = ttxfont_save;
	    }
	}
    }

/* Restore the correct color, if necessary */
    if (cur_color != cur_color_save)
    {
	cur_color = cur_color_save;
	need_devcolor = YES;
    }

/* Restore overlay mode */
    overlay = overlay_save;

/*
 * If they jump back into text they can continue right where they left off
 * (As long as they do nothing to move the pen between now and then.)
 */
    mov (last_widthr + font[ttxfont].dim[LETTER] * SIZE_FACTOR (path),
	 ytxshift);
    xold = ROUND (xold_f);
    yold = ROUND (yold_f);
    cfree ((char *) istring);
}

mov (hadd, vadd)
    double          hadd, vadd;
{
    xold_f += hadd * path_orient_dx + vadd * up_orient_dx;
    yold_f += hadd * path_orient_dy + vadd * up_orient_dy;
}

#define FONTCHECK_STRING "Vplot Binary fonT  \n"

load_font (ifont)
    int             ifont;
{
int             fd, length;
char            filename[120];
char            string[80];
char           *newfont;
int             offs[7];
static int      done[NUMGENFONT];
char           *stringptr;
int             fontcheck;
int             need_swab, file_ok;
int             nread;
int             vplotfontdirset;

    if (done[ttxfont])
    {
/*
 * Fonts are modulo NUMGENFONT. The first try in each slot determines
 * what happens. If the first try bombed, well, too bad, you get the
 * error font for all fonts that land there modulo NUMGENFONT forever more.
 */
	ttxfont = ERRFONT;
	return;
    }

/*
 * One way or another we're going to be done looking for this font
 * when we exit this routine.
 */
    done[ttxfont] = YES;


/*
 * Here is where we look for the font.
 * If the font is one of the "hardwired" ones, we look in the
 * SYSTEM_FONT_DIRECTORY for it by name [unless the environmental variable
 * "VPLOTFONTDIR" has been set in which case we look there instead].
 * If the font is not one of the "hardwired" ones, we look for
 * font number XX in the file "fontXX" in the current directory.
 *
 * In either case we also getpar for "fontXX" to see if the user has specified
 * an overriding file name on the command line.
 */
    vplotfontdirset = -1;

    if (ttxfont < NUM_FONTS)
    {
	if ((stringptr = getenv ("VPLOTFONTDIR")) != NULL)
	{
	    vplotfontdirset = YES;
	    sprintf (filename, "%s%s.bin", stringptr, font[ttxfont].name);
	}
	else
	{
	    vplotfontdirset = YES;
	    sprintf (filename, "%s%s.bin", VPLOT_FONT_DIR, font[ttxfont].name);
/*
	    vplotfontdirset = NO;
	    sprintf (filename, "%s%s.bin", SYSTEM_FONT_DIRECTORY, font[ttxfont].name);
*/
	}
    }
    else
    {
/*
 * The default place to look for a user-specified vplot binary font.
 */
	sprintf (filename, "./font%d.bin", ifont);
    }

/*
 * getpar parameter fontXX=name overrides all else. (Unless the font was
 * loaded at compile time in which case you don't get the opportunity.)
 */
    sprintf (string, "font%d", ifont);
    if (getpar (string, "s", filename))
	vplotfontdirset = -1;

/*
 * Try opening it.
 */
    if ((fd = open (filename, O_RDONLY)) == -1)
    {
/* Huh, it wasn't there. Complain and give up. */
	ERR (WARN, name,
	     "(gentext) Couldn't find font %d, file \"%s\".",
	     ttxfont, filename);
	if (vplotfontdirset == NO)
	{
	    ERR (COMMENT, name,
		 "(gentext) Perhaps you need to setenv VPLOTFONTDIR on this machine?");
	}
	else
	if (vplotfontdirset == YES)
	{
	    ERR (COMMENT, name,
		 "(gentext) Is your environmental variable VPLOTFONTDIR correct?");
	}
	else
	{
	    ERR (COMMENT, name,
	    "(gentext) Is the command-line parameter font%d=\"%s\" correct?",
		 ifont, filename);
	}
	ttxfont = ERRFONT;
	return;
    }


/*
 * At least there is a file there. But is it the right kind of file?
 * First check to make sure it has the magic sequence "Vplot Binary fonT  \n"
 * at the beginning. If it doesn't, it's junk, so don't read it in.
 */
    nread = read (fd, string, strlen (FONTCHECK_STRING));

    if (nread != strlen (FONTCHECK_STRING) ||
	strncmp (FONTCHECK_STRING, string, strlen (FONTCHECK_STRING)) != 0)
    {
	close (fd);
	ERR (WARN, name,
	     "(gentext) Font %d file \"%s\" is not a vplot font.",
	     ttxfont, filename);
	ttxfont = ERRFONT;
	return;
    }


/*
 * It begins promisingly... but is the opening string followed by
 * the binary integer FONTCHECK? Perhaps the file needs to be
 * byte-swapped?
 */
    file_ok = YES;
    need_swab = NO;

    nread = read (fd, (char *) &fontcheck, sizeof (int));

    if (nread != sizeof (int))
	file_ok = NO;
    else
    {
	if (fontcheck != FONTCHECK)
	{
	    swab_font ((char *) &fontcheck, (int) sizeof (int));
	    if (fontcheck != FONTCHECK)
		file_ok = NO;
	    else
	    {
		need_swab = YES;
#ifdef DEBUG
		ERR (WARN, name,
		     "(gentext) Font %d file \"%s\" needs byte-swapping.",
		     ttxfont, filename);
#endif
	    }
	}
    }

    if (file_ok == NO)
    {
font_garbled:
/*
 * I guess if I were really tidy I would free up the memory allocated
 * to a font found to be truncated on disk. This way the partial
 * font will remain in memory but marked as "unloaded" since it was
 * found to be partially damaged. Well, I guess this way a manic
 * debugger can examine it in memory...
 */
	close (fd);
	ERR (WARN, name,
	     "(gentext) Font %d file \"%s\" is garbled or truncated.",
	     ttxfont, filename);
	ttxfont = ERRFONT;
	return;
    }

/*
 *  Suck the fonts into memory, byte-swapping as we go if indicated.
 *  The start of the file contains the length and then the
 *  offsets from the beginning to the 7 structures defining the font.
 */
    nread = read (fd, (char *) &length, sizeof (int));
    if (nread != sizeof (int))
	goto font_garbled;
    if (need_swab)
	swab_font ((char *) &length, (int) sizeof (int));

    newfont = (char *) malloc ((unsigned) length);
    if (newfont == NULL)
    {
	close (fd);
	ERR (WARN, name,
	     "(gentext) Can't allocate memory for Font %d file \"%s\".",
	     ttxfont, filename);
	ttxfont = ERRFONT;
	return;
    }

/* Offsets to the 7 structures defining the font... */
    nread = read (fd, (char *) offs, 7 * sizeof (int));
    if (nread != 7 * sizeof (int))
	goto font_garbled;
    if (need_swab)
	swab_font ((char *) offs, 7 * (int) sizeof (int));

/* The font data itself */
    nread = read (fd, (char *) newfont, length);
    if (nread != length)
	goto font_garbled;
    if (need_swab)
	swab_font (newfont, length);

/* Done! */
    close (fd);

/* The vital parameters (dimensions, bounds, and lengths) */
    font[ttxfont].dim = (int *) (newfont + offs[0]);
/* Pointers to the addresses of the glyphs themselves */
    font[ttxfont].saddr = (int *) (newfont + offs[1]);
/* Left widths */
    font[ttxfont].swidthl = (int *) (newfont + offs[2]);
/* Right widths */
    font[ttxfont].swidthr = (int *) (newfont + offs[3]);
/* Vertical symbol hot-spot position */
    font[ttxfont].symbol = (int *) (newfont + offs[4]);
/* The actual glyph strokes */
    font[ttxfont].svec = (unsigned int *) (newfont + offs[5]);
/* Ligature data */
    font[ttxfont].lig = (int *) (newfont + offs[6]);

/* Whether or not this font has been loaded into memory or not yet */
    font[ttxfont].load = YES;
}

/*
 *
 *  source file:   ./filters/genlib/genvector.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 * Joe Dellinger (SOEST), June 23 1992
 *	Check the status of "lost" after every call to dev.plot,
 *	it might happen to change right after the first call of a
 *	vector pair.
 */

/*
 * Generic vector routine for devices that don't support fatness.
 * This version tries to be smart and minimize the "motion of the pen",
 * and tries to prolong strings of "draws" where possible.
 */
#include <stdio.h>
#include "../include/extern.h"
#define MOVE 0
#define DRAW 1

extern int      smart_clip;
extern int      lost;
extern int      fatvec ();

genvector (x1, y1, x2, y2, nfat, dashon)
    int             x1, y1, x2, y2;
    int             nfat, dashon;
{
static int      xlst, ylst;
int             d1, d2;

    if (nfat < 0)
	return;

    if (dashon)
    {
	dashvec (x1, y1, x2, y2, nfat, dashon);
	return;
    }

    if (nfat)			/* Recursively calls itself to make fat lines */
    {
	fatvec (x1, y1, x2, y2, nfat, dashon);
	return;
    }

    /*
     * Do clipping
     */
    if (!smart_clip)
	if (clip (&x1, &y1, &x2, &y2))
	    return;
    /*
     * Important special case: Zero-length vector at the end of what you've
     * already plotted. Don't need to do anything.
     */
    if (x1 == x2 && y1 == y2 && !lost && x1 == xlst && y1 == ylst)
    {
	return;
    }

    /*
     * Minimize movement of "pen"
     */
    if (!lost)
    {
	d1 = abs (x1 - xlst) + abs (y1 - ylst);
	d2 = abs (x2 - xlst) + abs (y2 - ylst);
	if (d2 < d1)
	{
	    d1 = x1;
	    d2 = y1;
	    x1 = x2;
	    y1 = y2;
	    x2 = d1;
	    y2 = d2;
	}
    }

    if ((x1 != xlst) || (y1 != ylst) || lost)
    {
	/* Make sure it is a move, not a draw */
	if (!lost && abs (x1 - xlst) <= 1 && abs (y1 - ylst) <= 1)
	{
	    /*
	     * We're within one pixel, so go ahead and draw a vector to the
	     * new point. This avoids having to leave and re-enter vector
	     * mode.
	     */
	    dev.plot (x1, y1, DRAW);

	    /*
	     * It is remotely possible that JUST NOW, when we did that DRAW,
	     * we caused the device to become "lost"... So check for that
	     * unfortunate unlikely case, and fix it. (I've since learned why
	     * a plethora of externals isn't good programming style...)
	     */
	    if (lost)
		dev.plot (x1, y1, MOVE);
	    /*
	     * If the device is STILL lost it's a BUG!
	     */
	}
	else
	{
	    dev.plot (x1, y1, MOVE);
	}
    }

    dev.plot (x2, y2, DRAW);
    xlst = x2;
    ylst = y2;
}

#endif
