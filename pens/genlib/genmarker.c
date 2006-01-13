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
	sprintf (txbuffer, "\\v%d ", type);
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
}

static void text_marker (char *txbuffer, int size, int npts, int *pvec)
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
