/* Bi-linear interpolation in 2-D */
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
/*^*/

#include "lint2.h"

static int    nx,  ny;
static float  ox,  oy;
static float  dx,  dy;
static float *xx, *yy;

void lint2_init( int n1, float o1, float d1 /* first grid axis */, 
		 int n2, float o2, float d2 /* second grid axis */, 
		 float* x, float* y         /* coordinates */) 
/*< initialize >*/
{
    nx = n1;
    ox = o1;
    dx = d1;
    xx = x;

    ny = n2;
    oy = o2;
    dy = d2;
    yy = y;
}

void lint2_lop( bool adj, bool add, int nm, int nd, float* mm, float*dd) 
/*< linear operator >*/
{
    int   ix, iy, im, id;
    float fx, fy;
    float gx, gy;
    float f;

    sf_adjnull (adj, add, nm, nd, mm, dd);
    
    for (id= 0; id< nd; id++) {

	f = (xx[id]-ox)/dx;
	ix = (int) f;  
	if ( ix < 0 || ix > nx-2) continue;
	fx=f-ix;
	gx= 1.-fx;
	
	f = (yy[id]-oy)/dy;
	iy = (int) f;  
	if ( iy < 0 || iy > ny-2) continue;
	fy=f-iy;
	gy= 1.-fy;
	
	im = nx*iy+ix;
	
	if (adj) {
	    mm[im  ]    += gx * gy * dd[id];
	    mm[im+1]    += fx * gy * dd[id];
	    mm[im+nx]   += gx * fy * dd[id];
	    mm[im+nx+1] += fx * fy * dd[id];
	} else {
	    dd[id] += 
		gx * gy * mm[im] + 
		fx * gy * mm[im+1] +
		gx * fy * mm[im+nx] +
		fx * fy * mm[im+nx+1];
	}

    }
}

