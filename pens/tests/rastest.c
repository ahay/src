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
 *  source file:   ./filters/Tests/rastest.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>

#include <rsf.h>
#include <rsfplot.h>

#define X 200
#define Y 200

int main (void)
{
    int             offset, xpix, ypix, bit, blast;
    float           xll, yll, xur, yur, ppi;
    unsigned char   **array;
    int             ii, jj;

    array = sf_ucharalloc2(X,Y);

    vp_init();
/*
 * Create a moire pattern
 */
    for (ii = 0; ii < Y; ii++)
    {
	for (jj = 0; jj < X; jj++)
	{
	    array[ii][jj] = ((((ii - 105) * (ii - 95) + 
			       (jj - 110) * (jj - 90)) / 77) % 14) + 1;
	    if (array[ii][jj] > 0)
	    {
		if (array[ii][jj] % 2 == 1)
		    array[ii][jj] = 0;
		array[ii][jj] /= 2;
	    }
	    else
	    {
/*
		if (array[ii][jj] % 2 == -1)
		    array[ii][jj] = 0;
*/
		array[ii][jj] /= -13;
	    }
	}
    }

    offset = 0;
    xpix = X;
    ypix = Y;
    bit = 0;
    xll = 0.;
    yll = 0.;
    xur = 10.;
    yur = 10.;
    ppi = 0;
    blast = 0;

/*
 * Decide what "style" plot this is to be.
 */
    vp_style (VP_STANDARD);

/*
 * Draw the raster.
 */
    vp_raster (array, bit, offset, 
	       xpix, ypix, xll, yll, xur, yur, 1);

/*
 * Draw a thin blue border around the whole thing.
 */
    vp_fat (0);
    vp_color (VP_BLUE);

/*
 * (xll,yll) is the lower-leftmost pixel of the raster plot,
 * and (xur-xll) is the width in vplot units and (yur-yll) is the
 * height.
 * If you think about this carefully, you'll see that then (xur,yur)
 * is not quite the upper-rightmost pixel of the raster plot, but is
 * one off.
 * Things were done this way because the height and width are
 * the real parameters that you want to specify.
 * Unfortunately it also means that the border isn't exactly symmetrical.
 */
    vp_move (xll, yll);
    vp_draw (xll, yur);
    vp_draw (xur, yur);
    vp_draw (xur, yll);
    vp_draw (xll, yll);

    exit(0);
}
