/* Simple tomographic operator */
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

#include <math.h>

#include <rsf.h>

#include "tomo.h"

/* velocity grid dimensions */
static int nz, nx, np, nq, dq;
static float dz, dx, dp, p0;

void tomo_init(int np1            /* number of slopes */,
	       int ns1, int ds1   /* depth steps */,
	       int n1, int n2     /* grid dimensions */, 
	       float d1, float d2 /* grid sampling */)
/*< initialize >*/
{
  nz = n1;
  nx = n2;
  dz = d1;
  dx = d2;

  nq = ns1;
  dq = ds1;

  np = np1;
  dp = 2./(np-1);
  p0 = -1.;
}

void tomo_lop(bool adj, bool add, int ns, int nt, float* s, float* t)
/*< linear operator >*/
{
    int is, ix, iz, ip, iy, it, i, kz;
    float p, x, deltax, distance;
    
    if (ns != nx*nz || nt != nx*np*nq) sf_error("%s: wrong size",__FILE__);
    
    /* initialize for lop
       if !add && adj, zero s
       if !add && !adj, zero t
       if add, do nothing */
    sf_adjnull(adj,add,ns,nt,s,t);
 
    for (i=0; i < nq; i++) { /* loop over initial depths */
		kz = i*dq;
		if (kz > nz-1) kz=nz-1;

		for (iy=0; iy < nx; iy++) { /* loop over sources */
			for (ip=0; ip < np; ip++) { /* loop over initial directions */
				p = p0 + ip*dp; /* initial slope */
				x = iy*dx; /* initial position (origin iz zero) */
	    
				is = iy*nz + kz; /* where in the grid is [kz,iy] */
				it = (i*nx + iy)*np + ip;
	    
				deltax = dz*p; /* shift in x */
				distance = hypotf(dz,deltax); /* Pythagor */
		
				for (iz=kz; iz > 0; iz--) { /* loop up in depth */
					if (adj) {
						s[is] += t[it]*distance;
					} else {
						t[it] += s[is]*distance; /* time is slowness * distance */
					}
		
					x += deltax;
					ix = 0.5 + x/dx; /* nearest point on the grid */
		    
					if (ix < 0 || ix >= nx) break; /* outside the grid */
		
					is = ix*nz + iz; 
				} /* iz */      
			} /* ip */ 
		} /* is */
    } /* i */
}
