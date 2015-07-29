/* 1-D internal convolution, adjoint is the input */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
 
#include <rsf.h>

#include "icai1.h"

static int nb, lg;
static float* bb;

void icai1_init (int na    /* filter length */, 
		 float* aa /* filter [na] */, 
		 int lag   /* filter lag (lag=1 is causal) */) 
/*< initialize >*/
{
    nb = na;
    bb = aa;
    lg = lag;
}

void icai1_lop (bool adj, bool add, int nx, int ny, float* xx, float* yy) 
/*< linear operator >*/
{
    int b, x, y;

    if(ny != nx) sf_error("%s: size problem: %d != %d",__FILE__,ny,nx);
    sf_adjnull (adj, add, nx, ny, xx, yy);
    
    for( b=0; b < nb; b++) {
	for( y = nb - lg; y <= ny - lg; y++) {
	    x = y - b + lg - 1;
	    if( adj) xx[x] += yy[y] * bb[b];
	    else     yy[y] += xx[x] * bb[b];
	}
    }
}

/* 	$Id: icai1.c 7107 2011-04-10 02:04:14Z ivlad $	 */
