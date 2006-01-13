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
 *  source file:   ./filters/utilities/dupside.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include "../include/vertex.h"
/*^*/

int dupside (register struct vertex *base)
/*< Determine if other sides in the polygon are
 * identical to the side specified by the
 * vertices v and v->b. >*/
{
register struct vertex *v;
register int    x1, x2, y1, y2;

    x1 = base->x;
    x2 = base->last->x;
    y1 = base->y;
    y2 = base->last->y;
    v = base->next;
    do
    {
	if (x1 == v->x && y1 == v->y && x2 == v->last->x && y2 == v->last->y)
	    return (1);
	if (x2 == v->x && y2 == v->y && x1 == v->last->x && y1 == v->last->y)
	    return (1);
	v = v->next;
    }
    while (v != base);
    return (0);
}
