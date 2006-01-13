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
 *  source file:   ./filters/genlib/nulldev.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

#include <stdio.h>

#include "../include/vertex.h"

void nulldev (void)
/*<  All purpose Do-nothing generic subroutine >*/
{ }

void nullclose (int status)
/*<  All purpose Do-nothing generic subroutine >*/
{ }

void nullarea (int npts, struct vertex  *head)
/*<  All purpose Do-nothing generic subroutine >*/
{ }

void nullraster(int xpix, int ypix, int xmin, int ymin, int xmax, int ymax, 
		unsigned char **raster_block, int orient, int dither_it)
/*<  All purpose Do-nothing generic subroutine >*/
{ }

void nullattributes(int command, int value, int v1, int v2, int v3)
/*<  All purpose Do-nothing generic subroutine >*/
{ }

int nullgetpoint (int *x, int *y)
/*<  All purpose Do-nothing generic subroutine >*/
{ return 0; }

int nullinteract (int what, FILE *controltty, char *string)
/*<  All purpose Do-nothing generic subroutine >*/
{ return 0; }

void nullvector(int x1, int y1, int x2, int y2, int nfat, int vpdashon)
/*<  All purpose Do-nothing generic subroutine >*/
{ }
