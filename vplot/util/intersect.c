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
 *  source file:   ./filters/utilities/intersect.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include	"../include/vertex.h"


intersect (x, crosses, head, scany)
    register int    x;
    register int   *crosses;
    struct vertex  *head;
    int             scany;
{
register struct vertex *v;
int             ncross, y;

    ncross = 0;
    v = head;
    if (scany)
    {
	do
	{
	    if (v->y > x && v->last->y > x)
		continue;
	    if (v->y < x && v->last->y < x)
		continue;

	    y = solve (x, v->y, v->x, v->last->y, v->last->x);
	    crosses[ncross++] = y;
	} while ((v = v->next) != head);
    }
    else
    {
	do
	{
	    if (v->x > x && v->last->x > x)
		continue;
	    if (v->x < x && v->last->x < x)
		continue;

	    y = solve (x, v->x, v->y, v->last->x, v->last->y);
	    crosses[ncross++] = y;
	} while ((v = v->next) != head);
    }
    return (ncross);
}
