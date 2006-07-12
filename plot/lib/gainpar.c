/* Finding clip and gpow. */
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

#include <math.h>

#include <rsf.h>
/*^*/

#include "gainpar.h"

static void gain (sf_file in, float **data, int n1, int n2,int step,
		  float pclip,float phalf,
		  float *clip, float *gpow, float bias, int nt, float* buf);

void vp_gainpar (sf_file in, float **data, 
		 int n1, int n2 /* panel size */,
		 int step       /* step on the first axis */,
		 float pclip    /* percentage clip */,
		 float phalf    /* percentage for gpow */,
		 float *clip    /* output clip */, 
		 float *gpow    /* output gpow */, 
		 float bias, 
		 int n3         /* number of panels */ ,
		 int panel      /* gain type */)
/*< Find clip and gpow parameters >*/
{
    int nt, i3, nclip, nhalf;
    float *buf, *clipnp, *gpownp;

    nt = n1 / step;
    buf = sf_floatalloc(nt*n2);
  
    if(panel >= 0) { /* gain from a particular panel */
	gain (in, data, n1, n2, step, pclip, phalf,
	      clip, gpow, bias, nt, buf);
    } else { /* gainpanel=each */
	clipnp = sf_floatalloc(n3);
	gpownp = sf_floatalloc(n3);

	for (i3=0; i3<n3; i3++) {
	    clipnp[i3] = 0.;
	    gpownp[i3] = 0.;	    
	    gain (in, data, n1, n2, step, pclip, phalf,
		  clipnp+i3, gpownp+i3, bias, nt, buf);
	    if (NULL == in) data += n1*n2; /* next panel */
	}

	nclip = SF_MAX(SF_MIN(n3*pclip/100. + .5,n3-1),0);
	nhalf = SF_MAX(SF_MIN(n3*phalf/100. + .5,n3-1),0);

	if (*clip==0.) *clip = sf_quantile (nclip,n3,clipnp);
	if (*gpow==0.) *gpow = sf_quantile (nhalf,n3,gpownp);
        free(clipnp);
        free(gpownp);
    }	

    free (buf);
}
 
static void gain (sf_file in, float **data, int n1, int n2, int step,
		  float pclip,float phalf, float *clipp, float *gpowp, 
		  float bias, int nt, float *buf)
{
    int j, i2, it, ntest, n;
    float clip, gpow, half;

    n = n2*nt;

    if (NULL != in) sf_floatread (data[0],n1*n2,in);

    for (j=i2=0; i2<n2; i2++) {
	for (it=0; it<nt; it++, j++) {
	    buf[j] = fabsf(data[i2][it*step]- bias);
	}
    }
    
    /* determine gain factors */
    clip = *clipp;

    if (clip==0.) {
	ntest = n*pclip/100. + .5;
	if (ntest >= n) { /* 100% clip */
	    clip = buf[0];		
	    for (j=1; j < n; j++) {
		if(buf[j] > clip) clip=buf[j];
	    }	
	} else {
	    clip = sf_quantile(SF_MAX(ntest,0),n,buf);
	}
	*clipp = clip;
    }
    
    if (*gpowp==0.) {
	ntest = n*phalf/100. + .5;
	if (ntest >= n) { /* 100% half */
	    half = buf[0];
	    for (j=1; j < n; j++) {
		if(buf[j] > half) half=buf[j];
	    }	
	} else {
	    half = sf_quantile (SF_MAX(ntest,0),n,buf);
	}
	
	if (clip==0. || half == clip || half/clip < 0.001) {
	    gpow = 1.;
	} else {
	    gpow = logf (.5) / logf (half / clip);
	    if (gpow < 0.1) {
		gpow = 0.1;
	    } else if (gpow > 10.) {
		gpow = 10.;
	    }
	}
	*gpowp = gpow;
    }
}
