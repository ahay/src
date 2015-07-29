/* Boundary masks for two component estimation */
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

#include "mask4freq.h"

void mask4freq (int nw         /* filter size */, 
	        int nj         /* dealiasing stretch */, 
	        int nx, int ny /* data size */, 
	        float *yy     /* data [ny][nx] */, 
	        bool *mm      /* mask [ny][nx]*/) 
/*< mask in 2-D  frequency filter>*/
{
    int ix, iy, iw, is, n, i;
    bool *xx;

    n=nx*ny;

    xx = sf_boolalloc(n);
    
    for (i=0; i < n; i++) {
	xx[i] = (bool) (yy[i] == 0.);
	mm[i] = false;
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj; ix < nx-nw*nj; ix++) {
	    i=ix + nx*iy;
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj;
		mm[i] = (bool) (mm[i] || xx[i-is]);
	    }
	}
    }
    
    free(xx);
}

/* 	$Id: mask4freq.c 5595 2010-03-21 16:54:14Z sfomel $	 */
