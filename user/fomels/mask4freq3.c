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

#include "mask4freq3.h"

void mask4freq3 (int nw         /* filter size */, 
	        int nj         /* dealiasing stretch */, 
	        int nx, int ny, int nz /* data size */, 
	        float *yy     /* data [ny][nx][nz] */, 
	        bool *mm      /* mask [ny][nx][nz] */) 
/*< mask in 2-D  frequency filter>*/
{
    int ix, iy,iz, iw, is, n, i;
    bool *xx;

    n=nx*ny*nz;

    xx = sf_boolalloc(n);
    
    for (i=0; i < n; i++) {
	xx[i] = (bool) (yy[i] == 0.);
	mm[i] = false;
    }
    
    for (iz=0; iz<nz; iz++) {
    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj; ix < nx-nw*nj; ix++) {
	    i=ix + nx*(iy+ny*iz);
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj;
		mm[i] = (bool) (mm[i] || xx[i-is]);
	    }
	}
    }
    }
    
    free(xx);
}

/* 	$Id$	 */
