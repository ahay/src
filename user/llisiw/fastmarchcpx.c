/* fast marching interface for marching from the surface. */
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

#include "fastmarchcpx.h"

static float *d;
static int *n, *in;

void fastmarchcpx_init(int *n_in    /* length */, 
		       float *d_in  /* sampling */) 
/*< initailize >*/
{
    int maxband;

    n = n_in;
    d = d_in;

    maxband = 0;
    if (n[0] > 1) maxband += 2*n[1]*n[2];
    if (n[1] > 1) maxband += 2*n[0]*n[2];
    if (n[2] > 1) maxband += 2*n[0]*n[1];

    sf_pqueue_init (10*maxband);

    in = sf_intalloc(n[0]*n[1]*n[2]);
}

void fastmarchcpx (float* time  /* time */,
		   float* t0    /* fixed traveltime */,
		   bool* m      /* known mask */,
		   float* v     /* slowness squared */)
/*< Run fast marching eikonal solver >*/
{
    float *p;
    int npoints, i;

    sf_pqueue_start();
    sf_neighbors_init (in, d, n, 1, time);

    /* initialize from boundary */
    for (npoints =  sf_neighbors_mask (v, t0, m ,true);
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
}    

void fastmarchcpx_close(void)
/*< free allocated storage >*/
{
    sf_pqueue_close();
    free(in);
}
