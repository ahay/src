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

#include <stdio.h>
#include <math.h>

#include <rsfplot.h>

#include "../include/params.h"
#include "../include/enum.h"
#include "../include/extern.h"
#include "../include/round.h"

#include "gentext.h"

static void text_marker (char *txbuffer, int size, int npts, int *pvec);

void genmarker (int npts, int type, int size, int *pvec)
/*< device-independent marker >*/
{
    int             savetxfont, savetxprec, savetxovly, savefat;
    struct s_txalign  savealign;
    extern float    fatmult_orig;
    char            txbuf[12];
    char           *txbuffer = txbuf;

    savealign.hor = txalign.hor;
    savealign.ver = txalign.ver;

    savetxfont = dev.txfont;
    savetxprec = dev.txprec;
    savetxovly = dev.txovly;
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
    dev.txovly = OVLY_NORMAL;

    if ((type > 32) && (type < 127))	/* is it a printing character? */
    {
	*txbuffer = (char) (type);
	*(txbuffer + 1) = '\0';
	text_marker (txbuffer, size, npts, pvec);
    }
    else
	if (type >= 127)		/* special non-ASCII character */
	{
	    sprintf (txbuffer, "\\v%d ", (int) type);
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
		    dev.txfont = MATH;
		    *txbuffer = (char) 57;
		    *(txbuffer + 1) = '\0';
		    text_marker (txbuffer, size, npts, pvec);
		    break;
		case 3:		/* '*' */
		    dev.txfont = MATH;
		    *txbuffer = (char) 33;
		    *(txbuffer + 1) = '\0';
		    text_marker (txbuffer, size, npts, pvec);
		    break;
		case 4:		/* circle */
		    dev.txfont = MISC;
		    *txbuffer = (char) 105;
		    *(txbuffer + 1) = '\0';
		    text_marker (txbuffer, size, npts, pvec);
		    break;
		case 5:		/* 'X' */
		    dev.txfont = MATH;
		    *txbuffer = (char) 60;
		    *(txbuffer + 1) = '\0';
		    text_marker (txbuffer, size, npts, pvec);
		    break;
		case 20:		/* square */
		    dev.txfont = MISC;
		    *txbuffer = (char) 72;
		    *(txbuffer + 1) = '\0';
		    text_marker (txbuffer, size, npts, pvec);
		    break;
		case 21:		/* triangle */
		    dev.txfont = MISC;
		    *txbuffer = (char) 73;
		    *(txbuffer + 1) = '\0';
		    text_marker (txbuffer, size, npts, pvec);
		    break;
		case 22:		/* diamond */
		    dev.txfont = MISC;
		    *txbuffer = (char) 74;
		    *(txbuffer + 1) = '\0';
		    text_marker (txbuffer, size, npts, pvec);
		    break;
		case 23:		/* star */
		    dev.txfont = MISC;
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
    dev.txfont = savetxfont;
    dev.txprec = savetxprec;
    dev.txovly = savetxovly;
    fat = savefat;
}

static void text_marker (char *txbuffer, int size, int npts, int *pvec)
{

    while (npts--)
    {
	xold = *pvec;
	yold = *(pvec + 1);
	if (dev.txfont < NUMGENFONT)
	{
	    gentext (txbuffer,
		     /* Character path direction */
		     (float) size * dev.aspect_ratio,
		     (float) 0,
		     /* Character up vector direction */
		     (float) 0,
		     (float) size);
	}
	else
	{
	    dev.text (txbuffer,
		      /* Character path direction */
		      (float) size * dev.aspect_ratio,
		      (float) 0,
		      /* Character up vector direction */
		      (float) 0,
		      (float) size);
	}
	pvec += 2;
    }
}
