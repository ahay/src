/* Boundary masks for plane-wave destruction */
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

#include "mask6.h"

void mask32 (bool both              /* left and right predictions */,
	     int nw                 /* filter size */, 
	     int nj1, int nj2       /* dealiasing stretch */, 
	     int nx, int ny, int nz /* data size */, 
	     float *yy              /* data [nz*ny*nx] */, 
	     bool **m               /* dip mask [both? 4:2][nz*ny*nx] */) 
/*< two-dip masks in 3-D >*/
{
    int ix, iy, iz, iw, is, i, n;
    bool *xx;

    n = nx*ny*nz;

    xx = sf_boolalloc(n);

    for (i=0; i < n; i++) {
	xx[i] = (bool) (yy[i] == 0.);
	m[0][i] = false;
	m[1][i] = false;
	if (both) {
	    m[2][i] = false;
	    m[3][i] = false;
	}
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=0; iy < ny-1; iy++) {
	    for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
		i = ix + nx * (iy + ny * iz);

		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj1;		  
		    m[0][i] = (bool) (m[0][i] || xx[i-is] || xx[i+nx+is]);
		}
	    }
	}
    }

    for (iz=0; iz < nz-1; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
		i = ix + nx * (iy + ny * iz);

		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj2;		  
		    m[1][i] = (bool) (m[1][i] || xx[i-is] || xx[i+ny*nx+is]);
		}
	    }
	}
    }

    if (!both) {
	free(xx);
	return;
    }

    for (iz=0; iz < nz; iz++) {
	for (iy=1; iy < ny; iy++) {
	    for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
		i = ix + nx * (iy + ny * iz);

		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj1;		  
		    m[2][i] = (bool) (m[2][i] || xx[i-is] || xx[i-nx+is]);
		}
	    }
	}
    }

    for (iz=1; iz < nz; iz++) {
	for (iy=0; iy < ny; iy++) {
	    for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
		i = ix + nx * (iy + ny * iz);

		for (iw = 0; iw <= 2*nw; iw++) {
		    is = (iw-nw)*nj2;		  
		    m[3][i] = (bool) (m[3][i] || xx[i-is] || xx[i-ny*nx+is]);
		}
	    }
	}
    }
	
    free(xx); 
}

void mask3 (int nw         /* filter size */, 
	    int nj         /* dealiasing stretch */, 
	    int nx, int ny /* data size */, 
	    float **yy     /* data */, 
	    bool **mm      /* mask */) 
/*< one-dip mask in 2-D >*/
{
    int ix, iy, iw, is;
    bool **xx;

    xx = sf_boolalloc2(nx,ny);
    
    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    xx[iy][ix] = (bool) (yy[iy][ix] == 0.);
	    mm[iy][ix] = false;
	}
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj; ix < nx-nw*nj; ix++) {
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj;
		mm[iy][ix] = (bool) (mm[iy][ix] || xx[iy+1][ix+is] || xx[iy][ix-is]);
	    }
	}
    }
    
    free(xx[0]);
    free(xx);
}

void mask6 (int nw           /* filter size */, 
	    int nj1, int nj2 /* dealiasing stretch */, 
	    int nx, int ny   /* data size */, 
	    float *yy       /* data [ny][nx] */, 
	    bool *mm        /* mask [ny][nx] */) 
/*< two-dip mask in 2-D >*/
{
    int ix, iy, iw, is, n, i;
    bool *xx;

    n = nx*ny;

    xx = sf_boolalloc(n);
    
    for (i=0; i < n; i++) {
	mm[i] = (bool) (yy[i] == 0.);
	xx[i] = false;
    }

    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj1; ix < nx-nw*nj1; ix++) {
	    i = ix + nx*iy;

	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj1;
		xx[i] = (bool) (xx[i] || mm[i+nx+is] || mm[i-is]);
	    }
	}
    }
    
    for (i=0; i < n; i++) {
	mm[i] = false;
    }
    
    for (iy=0; iy < ny-1; iy++) {
	for (ix = nw*nj2; ix < nx-nw*nj2; ix++) {
	    i = ix + nx*iy;
	    
	    for (iw = 0; iw <= 2*nw; iw++) {
		is = (iw-nw)*nj2;
		mm[i] = (bool) (mm[i] || xx[i+nx+is] || xx[i-is]);
	    }
	}
    }
    
    free(xx);
}
