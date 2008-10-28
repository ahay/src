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

#include <rsf.h>
#include <math.h>

#include "bsteer.h"

#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a < b ? a : b)

static float *beam; /* stacked beam */ 

void bsteer(float ***data,       
	    int nt, int nx, int ny, 
	    int npx, int npy,  
	    float opx, float dpx, float opy, float dpy, 
	    float dt, float dx, float dy, float ox, float oy,
	    bool **live, bool mode, int nlive,
	    float pmax, int lwind, int n1out,
            float xref, float yref,
            float ***semb)
/*< beam steering >*/
{
    float x,y;                     /* receiver positions */
    float px,py;                   /* slopes in x and y */
    int istart,istop;              /* time integration start and end indices */
    int ix,iy,ipx,ipy;             /* counters for position and slope */
    int i1;                        /* time counter or time window counter */
    int j1;                        /* time counter inside time window */
    int itshift;                   /* time shift in samples */
    float tshift;                  /* time shift */
    float p;                       /* slope amplitude */

    beam = sf_floatalloc(nt);

    /* loop over 2-D slopes */
    for (ipy = 0; ipy < npy; ipy++) {
	py = opy + ipy*dpy;

	for (ipx = 0; ipx < npx; ipx++) {
	    px = opx + ipx*dpx;

	    for (i1 = 0; i1 < nt; i1++) {
		beam[i1] = 0.;
	    }
            for (i1 = 0; i1 < n1out; i1++) {
		semb[ipy][ipx][i1] = 0.;
	    }

            /* do not do this beam if p is too big */
	    p = sqrt(px*px + py*py);
	    if (p <= pmax) {

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
			    if ( (2*tshift) > ((2*itshift+1)*dt) ) itshift += 1;

			    istart = MAX(itshift,0);
			    istop = MIN(nt+itshift,nt);
                            /* loop over samples of beam */
			    for (i1 = istart; i1 < istop; i1++){
                                /* sum the shifted input into the beam */
				beam[i1] += data[iy][ix][i1-itshift];
			    }
			}
		    }
		}

                /* normalize stack */
		for (i1 = 0; i1 < n1out; i1++){
		    for (j1 = 0; j1 < lwind; j1++){
			semb[ipy][ipx][i1] += beam[i1*lwind+j1]/nlive;
		    }
		}

	    }
	}
    }

    free(beam);

}
 

