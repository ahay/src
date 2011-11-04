/* 3-D warping. */

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

#include "stretch4.h"
#include "warp2.h"

static int n1, n2, n3, nx, ny;
static map4 map1, map3 ;
static float *trace1, **trace2, ***str2, ***str3, ***slice1;

void warp3_init(int n1_in, float o1, float d1,
		int n2_in, float o2, float d2,
		int n3_in, float o3, float d3 /* output grid */,
		int nt, int ny_in, int nx_in  /* input  grid */,
		float eps                     /* regularization */)
/*< initialize >*/
{
    n1 = n1_in;
    n2 = n2_in;
    n3 = n3_in;

    ny = ny_in;
    nx = nx_in;

    map1 = stretch4_init (n1, o1, d1, nt, eps);
    map3 = stretch4_init (n3, o3, d3, nx, eps);

    trace1 = sf_floatalloc(n1);
    trace2 = sf_floatalloc2(n2,n3);

    warp2_init(n2,o2,d2,
	       n3,o3,d3,
	       ny,nx,eps);
    
    str2   = sf_floatalloc3(ny,nx,n1);
    str3   = sf_floatalloc3(ny,nx,n1);
    slice1 = sf_floatalloc3(ny,nx,n1);
}

void warp3_close(void)
/*< free allocated storage >*/
{
    stretch4_close(map1);

    free(trace1);

    free(*trace2);
    free(trace2);

    free(**str2);
    free(*str2);
    free(str2);

    free(**str3);
    free(*str3);
    free(str3);

    free(**slice1);
    free(*slice1);
    free(slice1);
}

void warp3(float ***slice  /* [nx][ny][nt] input */,
	   float ***coord1 /* [nx][ny][nt] coordinates */,
	   float ***coord2 /* [nx][ny][nt] coordinates */,
	   float ***coord3 /* [nx][ny][nt] coordinates */,
	   float ***slice2 /* [n3][n2][n1] output */)
/*< apply warping >*/
{
    int i1, i2, i3;

    for (i3=0; i3 < nx; i3++) {
	for (i2=0; i2 < ny; i2++) {
	    stretch4_define (map1,coord1[i3][i2]);	    
	
	    stretch4_apply  (map1,slice[i3][i2],trace1);	
	    for (i1=0; i1 < n1; i1++) {
		slice1[i1][i3][i2] = trace1[i1];
	    }
	
	    stretch4_apply  (map1,coord2[i3][i2],trace1);
	    for (i1=0; i1 < n1; i1++) {
		str2[i1][i3][i2] = trace1[i1];
	    }

	    stretch4_apply  (map1,coord3[i3][i2],trace1);
	    for (i1=0; i1 < n1; i1++) {
		str3[i1][i3][i2] = trace1[i1];
	    }
	}
    }
    
    for (i1=0; i1 < n1; i1++) {
	warp2(slice1[i1],str2[i1],str3[i1],trace2);

	for (i3=0; i3 < n3; i3++) {
	    for (i2=0; i2 < n2; i2++) {
		slice2[i3][i2][i1] = trace2[i3][i2];
	    }
	}
    }
}

void fwarp3(float ***slice2 /* [n3][n2][n1] input */,
	    float ***coord1 /* [nx][ny][nt] coordinates */,
	    float ***coord2 /* [nx][ny][nt] coordinates */,
	    float ***coord3 /* [nx][ny][nt] coordinates */,
	    float ***slice  /* [nx][ny][nt] output */)
/*< apply forward warping >*/
{
    int i1, i2, i3;

    for (i3=0; i3 < nx; i3++) {
	for (i2=0; i2 < ny; i2++) {
	    stretch4_define (map1,coord1[i3][i2]);	    
	
	    stretch4_apply  (map1,coord2[i3][i2],trace1);
	    for (i1=0; i1 < n1; i1++) {
		str2[i1][i3][i2] = trace1[i1];
	    }

	    stretch4_apply  (map1,coord3[i3][i2],trace1);
	    for (i1=0; i1 < n1; i1++) {
		str3[i1][i3][i2] = trace1[i1];
	    }
	}
    }

    for (i1=0; i1 < n1; i1++) {
	for (i3=0; i3 < nx; i3++) {
	    for (i2=0; i2 < ny; i2++) {
		trace2[i3][i2] = slice2[i3][i2][i1];
	    }
	}

	fwarp2(trace2,str2[i1],str3[i1],slice1[i1]);
    }
    
    for (i3=0; i3 < nx; i3++) {
	for (i2=0; i2 < ny; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		trace1[i1] = slice1[i1][i3][i2];
	    }
	    
	    stretch4_define (map1,coord1[i3][i2]);	    	
	    stretch4_invert  (map1,slice[i3][i2],trace1);
	}
    }
}
