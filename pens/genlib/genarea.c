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
 *  source file:   ./filters/genlib/genarea.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include "../include/vertex.h"
#include "../include/params.h"
#include "../include/extern.h"

#include "polysubs.h"
#include "polyfix.h"

extern int      smart_clip;

extern int      Allgone;

void genarea (int npts, struct vertex  *head)
/*< Device Independent Polygon treatment.
 * Do a first-pass sort of clipping using
 * polysubs, and then finish the job by calling polyfix and
 * polystart. >*/
{
struct vertex  *v;
int             firstpoint, i;

    Allgone = 1;		/* Assume none left unless polyfix tells us */
    if (!smart_clip)
    {
	firstpoint = 2;
	xminclip (0, 0, &firstpoint);	/* Tell them all to get ready */
    }

    firstpoint = 1;

    v = head;
    for (i = 0; i < npts; i++)
    {
	if (!smart_clip)
	{
	    xminclip (v->x, v->y, &firstpoint);
	}
	else
	{
	    polyfix (v->x, v->y, &firstpoint);
	}
	v++;
    }
    if (!smart_clip)
    {
	firstpoint = -1;	/* Means this was the last point! */
	xminclip (0, 0, &firstpoint);
    }
    if (Allgone == 0)		/* If still 1, means there's nothing left! */
    {
	polystart ();
    }
}
