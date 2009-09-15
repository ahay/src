/* 1D denoising using local polynomial fitting. */
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
#include "multidivn1.h"
#include "weight2.h"

int main (int argc, char* argv[]) 
{
    int n1, n2, nfw, n[SF_MAX_DIM], m[SF_MAX_DIM], rec[SF_MAX_DIM], i1, i2, k, i3, p, rect, niter;
    bool boundary, verb;
    
    float *trace, *tempt, *temp, *d, *fixd, *g, *f, *pre, mean, tindex, err;
    sf_filter aa;
    sf_file in, out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);
    /* get the trace length (n1) and the number of traces (n2) and n3*/
    
    if (!sf_getint("nfw",&nfw)) sf_error("Need integer input");
    /* filter-window length (positive and odd integer) */

    if (!sf_getint("rect",&rect)) sf_error("Need integer input");
    /* local smoothing radius */

    if (!sf_getbool("boundary",&boundary)) boundary = false;
    /* if y, boundary is data, whereas zero*/

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    
    trace = sf_floatalloc(n1);
    tempt = sf_floatalloc(n1);
    temp = sf_floatalloc(nfw);
    fixd = sf_floatalloc(nfw*2);
    d = sf_floatalloc(nfw*2);
    g = sf_floatalloc(nfw);
    f = sf_floatalloc(nfw*2);
    pre = sf_floatalloc(nfw*2);

    for(i3=0; i3 < nfw; i3 ++) {
	fixd[0*nfw+i3] = 1.;
	fixd[1*nfw+i3] = (i3+1.)/sqrtf((nfw+1.)*(2.*nfw+1.)/6.); /* normalize by RMS=1 (sqrt(sum_N(t-tau)^2)/N=1) */
    }
 
    aa = NULL;

    for(i3=0; i3 < SF_MAX_DIM; i3 ++) {
	n[i3] = 1;
	m[i3] = 1;
	rec[i3] = 1;
    }
    n[0] = nfw;
    n[1] = 2;
    m[0] = nfw;
    rec[0] = rect;

    multidivn_init(2, 1, nfw, m, rec, d, aa, verb);

    mean = 0.;
    for (i3=0; i3 < nfw*2; i3++) {
	mean += fixd[i3]*fixd[i3];
    }
    mean = sqrtf (mean/(nfw*2));

    for(i3=0; i3< nfw*2; i3++) {
	fixd[i3] /= mean;
    }
    
    for(i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);
	
	for(i1=0; i1 < n1; i1++){
	    tempt[i1]=trace[i1];
	}
	
	/************1D local polynomail fitting****************/

	for(i1=0; i1 < n1; i1++) {
	    for(p=0; p < nfw*2; p++) {
		d[p] = fixd[p];
	    }	    
	    for(p=0; p < nfw; p++) {
		for(k=0; k < nfw; k++) {
		    if ((i1+k+p-nfw+1)<0) {
			if(boundary) {
			    g[k] = tempt[0];
			} else {
			    g[k] = 0.;
			}
		    } else if((i1+k+p-nfw+1)>(n1-1)) {
			if(boundary) {
			    g[k] = tempt[n1-1];
			} else {
			    g[k] = 0.;
			}
		    } else {
			g[k] = tempt[i1+k+p-nfw+1];
		    }
		}
		for(k=0; k < nfw; k++) {
		    temp[k] = g[k];
		}
		
		for (k=0; k < nfw; k++) {
		    g[k] /= mean;
		}
		multidivn (g,f,niter);
		for (k=0; k < nfw*2; k++) {
		    d[k] *= mean;
		}
		weight2_lop(false,false,nfw*2,nfw,f,g);
		pre[1*nfw+p] = g[nfw-1-p];
		err = 0.;
		for (k=0; k < nfw; k++) {
		    err += sqrtf((g[k]-temp[k])*(g[k]-temp[k]));
		}
		pre[0*nfw+p] = err;
		for (k=0; k < nfw*2; k++) {
		    d[k] = fixd[k];
		}		    
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
	    trace[i1] = pre[1*nfw+nfw-1];
	}

	sf_floatwrite(trace,n1,out);
    }
    
    exit (0);
}

/* 	$Id$	 */
