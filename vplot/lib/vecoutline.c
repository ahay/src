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
 *  source file:   ./filters/utilities/vecoutline.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>
#include "vertex.h"
#include "device.h"

void vecoutline (device dev, struct vertex  *head)
{
    /*
     * Draw the outline of the polygon pointed to by 'head'.  If any side of
     * the polygon is identical to any other (that is they have identical
     * endpoints), then neither line is drawn.  This allows among things a
     * doughnut to be defined by a single polygon without the line that
     * connects the inner and outer being plotted. 
     */
    int    xlast, ylast;
    struct vertex *v;

    xlast = head->last->x;
    ylast = head->last->y;
    v = head;
    do
    {
	if (!dupside (v))
	    dev->vector (v->x, v->y, xlast, ylast, dev->afat, 0);
	xlast = v->x;
	ylast = v->y;
    } while ((v = v->next) != head);
}
