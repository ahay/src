/* 1-D non-stationary smoothing as a linear operator. */
/*
  Copyright (C) 2016 University of Texas at Austin

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

static int n1, n2;
static float **rct, *trace, *sft;
static sf_ntriangle tr;

void nsmooth1_init(int ns /* number of samples in a trace */,
		   int nt /* number of traces */,
		   float **rect /* [n2][n1] smoothing radius */)
/*< initialize >*/
{
    int rect1, i1, i2;

    n1 = ns;
    n2 = nt;
    rct = rect;
    
    sft = sf_floatalloc(n1);

    rect1=1;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (rct[i2][i1] > rect1) rect1=ceilf(rct[i2][i1]);
	    sft[i1] = 0.0f;
	}
    }

    trace = sf_floatalloc(n1);
    tr = sf_ntriangle_init (rect1,n1);
}

void nsmooth1_close(void)
/*< clean allocated storage >*/
{
    free(sft);
    free(trace);
    sf_ntriangle_close(tr);
}

void nsmooth1_lop(bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i1, i2;

    if (nx != n1*n2 || ny != n1*n2) sf_error("%s: Wrong size",__FILE__);

    sf_adjnull(adj,add,nx,ny,x,y);
    
    for (i2=0; i2 < n2; i2++) {	
	if (adj) {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1]=y[i1+i2*n1];
	    }
	    sf_nsmooth2 (tr,0,1,false,rct[i2],sft,trace);
	    for (i1=0; i1 < n1; i1++) {
		x[i1+i2*n1] += trace[i1];
	    }
	} else {
	    for (i1=0; i1 < n1; i1++) {
		trace[i1]=x[i1+i2*n1];
	    }
	    sf_nsmooth (tr,0,1,false,rct[i2],sft,trace);
	    for (i1=0; i1 < n1; i1++) {
		y[i1+i2*n1] += trace[i1];
	    }
	}
    }
}
