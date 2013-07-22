/* Fast marching for reverse time. */
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

#include "timecont.h"

void timecont (float* time                /* time */, 
	       float* t0                  /* time at the surface */,
	       float* v                   /* slowness */, 
	       int* in                    /* in/front/out flag */, 
	       int   n3,  int n2,  int n1 /* dimensions */,
	       float d3,float d2,float d1 /* sampling */,
	       int order                  /* accuracy order (1,2,3) */,
	       bool forwd                 /* forward or backward */)
/*< Run fast marching eikonal solver >*/
{
    float d[3], *p;
    int n[3], npoints, i, maxband;

    maxband = 0;
    if (n1 > 1) maxband += 2*n2*n3;
    if (n2 > 1) maxband += 2*n1*n3;
    if (n3 > 1) maxband += 2*n1*n2;

    sf_pqueue_init (10*maxband);
    
    n[0] = n1; d[0] = d1;
    n[1] = n2; d[1] = d2;
    n[2] = n3; d[2] = d3;

    sf_pqueue_start();
    sf_neighbors_init (in, d, n, order, time);

    if (forwd) {
	for (npoints =  sf_neighbors_surface (v, t0, true);
	     npoints > 0;
	     npoints -= sf_neighbours(i)) {
	    /* Pick smallest value in the NarrowBand
	       mark as good, decrease points_left */
	
	    p = sf_pqueue_extract();

	    if (p == NULL) {
		sf_warning("%s: heap exausted!",__FILE__);
		break;
	    }
	
	    i = p - time;

	    in[i] = SF_IN;
	}
    } else {
	for (npoints =  sf_neighbors_surface (v, t0, false);
	     npoints > 0;
	     npoints -= sf_neighbours2(i)) {
	    /* Pick largest value in the NarrowBand
	       mark as good, decrease points_left */
	
	    p = sf_pqueue_extract2();

	    if (p == NULL) {
		sf_warning("%s: heap exausted!",__FILE__);
		break;
	    }
	
	    i = p - time;

	    in[i] = SF_IN;
	}
    }
    
    sf_pqueue_close();
}

/* 	$Id: fastmarch.c 1507 2005-10-22 04:01:28Z savap $	 */
