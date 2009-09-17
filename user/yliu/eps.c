/* 1D edge-preserving smoothing (EPS). */
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

#include <math.h>
#include "deviation.h"
#include "mean.h"

float eps(float *inp     /* input data */,
	  float *outp    /* output data */,
	  int n          /* data size */,
	  int nfw        /* filter size */,
	  bool boundary  /* boundary flag */)
/*< edge-preserving smoothing >*/
{
    float *tempt, *g, *pre, tindex;
 
    int i1, k, p;

    tempt = sf_floatalloc(n);
    g = sf_floatalloc(nfw);
    pre = sf_floatalloc(nfw*2);

    for(i1=0; i1 < n; i1++){
	tempt[i1] = inp[i1];
    }

    for(i1=0; i1 < n; i1++) {
	for(p=0; p < nfw; p++) {
	    for(k=0; k < nfw; k++) {
		if ((i1+k+p-nfw+1)<0) {
		    if(boundary) {
			g[k] = tempt[0];
		    } else {
			g[k] = 0.;
		    }
		} else if((i1+k+p-nfw+1)>(n-1)) {
		    if(boundary) {
			g[k] = tempt[n-1];
		    } else {
			g[k] = 0.;
		    }
		} else {
		    g[k] = tempt[i1+k+p-nfw+1];
		}
	    }
	    pre[1*nfw+p] = mean(g,nfw);
	    pre[0*nfw+p] = sdeviation(g,nfw);
	}
	for(p=0; p < (nfw-1); p++) {
	    if(pre[0*nfw+p]<pre[0*nfw+p+1]){
		tindex=pre[0*nfw+p];
		pre[0*nfw+p]=pre[0*nfw+p+1];
		pre[0*nfw+p+1]=tindex;
		tindex=pre[1*nfw+p];
		pre[1*nfw+p]=pre[1*nfw+p+1];
		pre[1*nfw+p+1]=tindex;
	    }
	}
	outp[i1] = pre[1*nfw+nfw-1];
    }
    return outp[(n-1)/2];
}

/* 	$Id$	 */
