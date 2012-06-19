/* Fast sweeping. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include "sweeping.h"

void sweep (int niter                  /* maximum number of iterations */,
	    float* time                /* time */, 
	    float* v                   /* slowness squared */, 
	    int   n3,  int n2,  int n1 /* dimensions */,
	    float o3,float o2,float o1 /* origin */,
	    float d3,float d2,float d1 /* sampling */,
	    float s3,float s2,float s1 /* source */,
	    int   b3,  int b2,  int b1 /* box around the source */)
/*< Run fast sweeping >*/
{
    int iter, sw, j1, j2, j3;
    

    for(iter=0; iter < niter; iter++) {
	sf_warning("iter=%d",iter);

	for (sw=0; sw < 8; sw++) {
	    j1 = (sw & 1)? -1:1;
	    j2 = (sw & 2)? -1:1;
	    j3 = (sw & 4)? -1:1;	    
	}
    }
}
