/* 1-D internal convolution, adjoint is the filter */
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

#include "icaf1.h"

static int nx, lag;
static float* xx;

void icaf1_init (int ny    /* filter length */, 
		 float* yy /* data [nx] */, 
		 int lag1  /* filter lag (lag=1 is causal) */) 
/*< initialize >*/
{
    nx = ny;
    xx = yy;
    lag = lag1;
}

void icaf1_lop (bool adj, bool add, int nb, int ny, float* bb, float* yy) 
/*< linear operator >*/
{
    int b, x, y;

    if(ny != nx) sf_error("%s: size problem: %d != %d",__FILE__,ny,nx);
    sf_adjnull (adj, add, nb, ny, bb, yy);
    
    for( b=0; b < nb; b++) {
	for( y = SF_MAX(lag,b+1); y <= ny; y++) { x = y - b - 1;
	    if( adj) bb[b] += yy[y-lag] * xx[x];
	    else     yy[y-lag] += bb[b] * xx[x];
	}
    }
}

/* 	$Id: icaf1.c 838 2004-10-25 11:10:38Z fomels $	 */
