/* 1D edge-preserving local polynomial fitting (ELPF). */
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

#include "elpf.h"
#include <math.h>

static float *fixd, *d, mean;
static int m[SF_MAX_DIM], rec[SF_MAX_DIM];

void elpf_init (int n          /* data size */,
		int nfw        /* filter size */,
		int rect       /* smoothing radius */,
		bool verb      /* verbosity flag */)
/*< initialize >*/
{
    int i, n2;
    sf_filter aa;

    n2 = nfw*2;

    fixd = sf_floatalloc(n2);
    d = sf_floatalloc(n2);

    for(i=0; i < nfw; i ++) {
	fixd[0*nfw+i] = 1.;
	fixd[1*nfw+i] = (i+1.)/sqrtf((nfw+1.)*(2.*nfw+1.)/6.); 
    } /* normalize by RMS=1 (sqrt(sum_N(t-tau)^2)/N=1) */
 
    aa = NULL;

    for(i=0; i < SF_MAX_DIM; i ++) {
	m[i] = 1;
	rec[i] = 1;
    }
    m[0] = nfw;
    rec[0] = rect;

    sf_multidivn_init(2, 1, nfw, m, rec, d, aa, verb);

    mean = 0.;
    for (i=0; i < n2; i++) {
	mean += fixd[i]*fixd[i];
    }
    mean = sqrtf (mean/n2);

    for(i=0; i< n2; i++) {
	fixd[i] /= mean;
    }
}

void elpf_close (void)
/*< free allocated storage >*/
{
    free (fixd);
    free (d);
    sf_multidivn_close();   
}

float elpf (float *inp     /* input data */,
	    float *outp    /* output data */,
	    int n          /* data size */,
	    int nfw        /* filter size */,
	    int niter      /* iteration number */,
	    bool boundary  /* boundary flag */)
/*< local polynomial fitting >*/
{
    float *tempt, *temp, *g, *f, *pre, tindex, err;
 
    int i1, k, p;

    tempt = sf_floatalloc(n);
    temp = sf_floatalloc(nfw);
    g = sf_floatalloc(nfw);
    f = sf_floatalloc(nfw*2);
    pre = sf_floatalloc(nfw*2);

    for(i1=0; i1 < n; i1++){
	tempt[i1] = inp[i1];
    }

    for(i1=0; i1 < n; i1++) {
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

	    for(k=0; k < nfw; k++) {
		temp[k] = g[k];
	    }
	    
	    for (k=0; k < nfw; k++) {
		g[k] /= mean;
	    }

	    sf_multidivn (g,f,niter);

	    for (k=0; k < nfw*2; k++) {
		d[k] *= mean;
	    }
	    sf_weight2_lop(false,false,nfw*2,nfw,f,g);
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
	outp[i1] = pre[1*nfw+nfw-1];

    }
    return outp[(n-1)/2];
}

/* 	$Id$	 */
