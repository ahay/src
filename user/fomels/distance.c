/* Distance computation by fast marching. */
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

#include "distance.h"

void distance_init (int n3,int n2,int n1, int np) 
/*< Initialize data dimensions >*/
{
    int maxband;
    
    maxband = 0;
    if (n1 > 1) maxband += 2*n2*n3;
    if (n2 > 1) maxband += 2*n1*n3;
    if (n3 > 1) maxband += 2*n1*n2;

    sf_pqueue_init (100*SF_MAX(maxband,np));
}

void distance (int np         /* number of points */, 
	       float **points /* point coordinates [np][3] */,
	       float* dist    /* distance */, 
	       float* v       /* slowness squared */,
	       int* in                    /* in/front/out flag */, 
	       int n3,int n2,int n1       /* dimensions */,
	       float o3,float o2,float o1 /* origin */,
	       float d3,float d2,float d1 /* sampling */,
	       int order                  /* accuracy order (1,2,3) */)
/*< Run fast marching eikonal solver >*/
{
    float d[3], o[3], *p;
    int n[3], npoints, i;
    
    n[0] = n1; o[0] = o1; d[0] = d1;
    n[1] = n2; o[1] = o2; d[1] = d2;
    n[2] = n3; o[2] = o3; d[2] = d3;

    sf_pqueue_start();
    sf_neighbors_init (in, d, n, order, dist);

    for (npoints =  sf_neighbors_distance (np, v, points, d, o);
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

	i = p - dist;

	in[i] = SF_IN;
    }
}

void distance_close (void)
/*< Free allocated storage >*/
{
    sf_pqueue_close();
}

/* 	$Id: distance.c 754 2004-08-24 07:59:16Z fomels $	 */
