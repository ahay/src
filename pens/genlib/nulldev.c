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

void nullmidpoly (int x, int y)
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

void nullplot(int x, int y, int draw)
/*<  All purpose Do-nothing generic subroutine >*/
{ }
