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
 *  source file:   ./filters/genlib/vecarea.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include <stdlib.h>

#include "../include/pat.h"
#include "../include/vertex.h"
#include "../include/params.h"
#include "../include/extern.h"

#include "../utilities/util.h"

void vecarea (int npts, struct vertex  *head)
/*< device-independent area >*/
{
register int    x, y, i;
register int    xstr, xend, ystr, yend;
int             skip;
int             ncross;
int             vminx, vmaxx, vminy, vmaxy;
struct vertex  *xhead, *yhead, *v;
int            *crosses;

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

    if (vmaxx > xwmax)
	vmaxx = xwmax;
    if (vminx < xwmin)
	vminx = xwmin;
    if (vmaxy > ywmax)
	vmaxy = ywmax;
    if (vminy < ywmin)
	vminy = ywmin;

    if ((pat[ipat] .ydim > 0) || (pat[ipat] .xdim == 1))
    {
	/* stretch polygon in y-direction */
	v = yhead;
	do
	{
	    v->y = 2 * (v->y) + 1;
	    v = v->next;
	} while (v != yhead);

	skip = (pat[ipat] .xdim == 1) ? 1 : pat[ipat] .ydim;
	for (y = vminy; y <= vmaxy; y += skip)
	{
	    ncross = intersect (2 * y, crosses, yhead, 1);
	    sort (crosses, ncross);
	    for (i = 0; i < ncross; i += 2)
	    {
		xstr = crosses[i];
		xend = crosses[i + 1];
		dev.vector (xstr, y, xend, y, 0, 0);
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

    if ((pat[ipat] .xdim > 1) && (pat[ipat] .ydim > 1))
    {
	/*
	 * expand in x 
	 */
	v = xhead;
	do
	{
	    v->x = 2 * v->x + 1;
	    v = v->next;
	} while (v != xhead);

	skip = pat[ipat] .xdim;
	for (x = vminx; x <= vmaxx; x += skip)
	{
	    ncross = intersect (2 * x, crosses, xhead, 0);
	    sort (crosses, ncross);
	    for (i = 0; i < ncross; i += 2)
	    {
		ystr = crosses[i];
		yend = crosses[i + 1];
		dev.vector (x, ystr, x, yend, 0, 0);
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
    }
    free ((char *) crosses);
}
