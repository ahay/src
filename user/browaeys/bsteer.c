/* Beam steering */
/*
  Copyright (C) 2008 University of Texas at Austin
  Adapted from Steve Cole, Stanford University, 1995
  
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
/*^*/

#include "bsteer.h"
/*^*/

#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a < b ? a : b)

void bsteer(float ***data,       
	    int nt, int nx, int ny, 
	    int npx, int npy,  
	    float opx, float dpx, float opy, float dpy, 
	    float dt, float dx, float dy, float ox, float oy,
	    bool **live, bool mode, int nlive,
            float xref, float yref,
            float **semb)
/*< beam steering >*/
{
    float x,y;                     /* receiver positions */
    float px,py;                   /* slopes in x and y */
    int istart,istop;              /* time integration start and end indices */
    int ix,iy,ipx,ipy;             /* counters for position and slope */
    int it;                        /* time counter */
    int itshift;                   /* time shift in samples */
    float tshift;                  /* time shift */
    float beam;                    /* beam stack */

    /* loop over 2-D slopes */
    for (ipy = 0; ipy < npy; ipy++) {
	py = opy + ipy*dpy;

	for (ipx = 0; ipx < npx; ipx++) {
	    px = opx + ipx*dpx;

            /* clear arrays */
	    semb[ipy][ipx] = 0.;
	    beam = 0.;
	    
	    /* loop over receivers */
	    for (iy = 0; iy < ny; iy++) {

		for (ix = 0; ix < nx; ix++) {

		    if (live[iy][ix]) {

                        /* position compared to reference trace */
			y = iy*dy + oy - yref;
			x = ix*dx + ox - xref;

                        /* compute the necessary time shift to align */
			/* the current trace with the reference trace */
			if (mode) {
			    tshift = px*x + py*y;
			} else {
                            /* py is apparent slowness, px is azimuth */
			    tshift = py*( cosf(SF_PI*px/180.)*x + sinf(SF_PI*px/180.)*y );
			}

                        /* Nearest integer */
			itshift = floorf(tshift/dt);
			if ( (2.*tshift) > ((2*itshift+1)*dt) ) itshift += 1;

                        /* loop over samples of beam */
			istart = MAX(itshift,0);
			istop = MIN(nt+itshift,nt);
			for (it = istart; it < istop; it++){
                            /* sum the shifted input into the beam */
			    beam += data[iy][ix][it-itshift];
			}

		    }
		}
	    }

            /* normalize stack */
	    semb[ipy][ipx] = beam/nlive;

	}
    }

}
 

