/* Fast marching main interface. */
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

#include "fastmarch.h"

static int *n, order, *in;
static float *o, *d;

void fastmarch_init (int* n1    /* dimensions */,
		     float* o1  /* origin */,
		     float* d1  /* sampling */,
		     int order1 /* accuracy order */) 
/*< Initialize data dimensions >*/
{
    int maxband;
    
    n = n1;
    order = order1;
    o = o1;
    d = d1;

    in = sf_intalloc(n[0]*n[1]*n[2]);

    maxband = 0;
    if (n[0] > 1) maxband += 2*n[1]*n[2];
    if (n[1] > 1) maxband += 2*n[0]*n[2];
    if (n[2] > 1) maxband += 2*n[0]*n[1];

    sf_pqueue_init (10*maxband);
}

void fastmarch (float* time                /* time */, 
		float* v                   /* slowness squared */, 
		float* s                   /* source */)
/*< Run fast marching eikonal solver >*/
{
    float xs[3], *p;
    int b[3], npoints, i;
    bool plane[3];

    xs[0] = s[0]-o[0]; b[0] = 1; plane[0] = false;
    xs[1] = s[1]-o[1]; b[1] = 1; plane[1] = false;
    xs[2] = s[2]-o[2]; b[2] = 1; plane[2] = false; 

    sf_pqueue_start();
    sf_neighbors_init (in, d, n, order, time);

    for (npoints =  sf_neighbors_nearsource (xs, b, d, v, plane);
	 npoints > 0;
	 npoints -= sf_neighbours(i)) {
	/* Pick smallest value in the NarrowBand
	   mark as good, decrease points_left */

	/* sf_warning("npoints=%d",npoints); */

	p = sf_pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p - time;

	in[i] = SF_IN;
    }
}

void fastmarch_close (void)
/*< Free allocated storage >*/
{
    sf_pqueue_close();
}

/* 	$Id: fastmarch.c 5686 2010-04-07 16:33:34Z llisiw $	 */
