/* Fast marching on spherical coordinates. */
/*
  Copyright (C) 2022 University of Texas at Austin
  
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

#include "fastmarchrtp.h"

void fastmarchrtp_init (int n3,int n2,int n1) 
/*< Initialize data dimensions >*/
{
    int maxband;
    
    maxband = 0;
    if (n1 > 1) maxband += 2*n2*n3;
    if (n2 > 1) maxband += 2*n1*n3;
    if (n3 > 1) maxband += 2*n1*n2;
	sf_warning("maxband=%d",5*maxband);
    sf_pqueue_init (50*maxband);
}

void fastmarchrtp (float* time                /* time */, 
		float* v                   /* slowness squared */, 
		int* in                    /* in/front/out flag */, 
		bool* plane                /* if plane source */,
		int   n3,  int n2,  int n1 /* dimensions */,
		float o3,float o2,float o1 /* origin */,
		float d3,float d2,float d1 /* sampling */,
		float s3,float s2,float s1 /* source */,
		int   b3,  int b2,  int b1 /* box around the source */,
		int order                  /* accuracy order (1,2,3) */)
/*< Run fast marching eikonal solver >*/
{
    float xs[3], d[3], *p;
    int n[3], b[3], npoints, i;
    
    n[0] = n1; xs[0] = 0; b[0] = b1; d[0] = d1;
    n[1] = n2; xs[1] = 0; b[1] = b2; d[1] = d2;
    n[2] = n3; xs[2] = 0; b[2] = b3; d[2] = d3;

    sf_pqueue_start();
    sf_neighbors_init (in, d, n, order, time);

    for (npoints =  sf_neighbors_nearsource_rtp (xs, b, d, v, plane);
	 npoints > 0;
	 npoints -= sf_neighbours(i)) {
	/* Pick smallest value in the NarrowBand
	   mark as good, decrease points_left */

 	/*sf_warning("npoints=%d",npoints);*/

	p = sf_pqueue_extract();

	if (p == NULL) {
	    sf_warning("%s: heap exausted!",__FILE__);
	    break;
	}
	
	i = p - time;

	in[i] = SF_IN;
    }	
	
}

void fastmarchrtp_close (void)
/*< Free allocated storage >*/
{
    sf_pqueue_close();
}


