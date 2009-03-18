/* 2-D Laplacian operator */
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

#include "laplacian.h"

static int type,n1,n2;
static float d1,d2, center, corner;

void laplacian_init(int type1          /* operator type */,
		    int nz, int nx     /* dimensions */,
		    float dz, float dx /* sampling */)
/*< initialize >*/
{
    float s1, s2;

    type = type1;
    n1 = nz;
    n2 = nx;
    d1 = 1./(dz*dz);
    d2 = 1./(dx*dx);

    switch(type) {
	case 0:
	    center = -2.0*(d1+d2);
	    break;
	case 1:
	    corner = (d1+d2)/12.0;
	    s1 = (5*d1-d2)/6.0;
	    s2 = (5*d2-d1)/6.0;
	    d1 = s1;
	    d2 = s2;
	    center = -2.0*(2.0*corner+d1+d2);
	    break;
	default:
	    sf_error("%s: Unknown Laplacian type",__FILE__);
    }
}

void laplacian(float **uin  /* [nx][nz] */, 
	       float **uout /* [nx][nz] */)
/*< apply >*/
{
    int i1, i2;

    for (i2=0; i2 < n2; i2++) {
	uout[i2][0]    = 0.;
	uout[i2][n1-1] = 0.;
    }
    for (i1=0; i1 < n1; i1++) {
	uout[0][i1]    = 0.;
	uout[n2-1][i1] = 0.;
    }
    
    switch(type) {
	case 0:
    	    for (i2=1; i2 < n2-1; i2++) {
		for (i1=1; i1 < n1-1; i1++) {
		    uout[i2][i1] = 
			d1*(uin[i2][i1-1]+uin[i2][i1+1]) +
			d2*(uin[i2-1][i1]+uin[i2+1][i1]) +
			center*uin[i2][i1];
		}
	    }
	    break;
	case 1:
	    for (i2=1; i2 < n2-1; i2++) {
		for (i1=1; i1 < n1-1; i1++) {
		    uout[i2][i1] = 
			d1*(uin[i2][i1-1]+uin[i2][i1+1]) +
			d2*(uin[i2-1][i1]+uin[i2+1][i1]) +
			corner*(uin[i2+1][i1-1]+uin[i2+1][i1+1]  +
				uin[i2-1][i1-1]+uin[i2-1][i1+1]) +
			center*uin[i2][i1];
		}
	    }
	    break;
	default:
	    sf_error("%s: Unknown Laplacian type",__FILE__);
    }
}

