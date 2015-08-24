/* bilateral filtering. */
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

#include <stdio.h>
#include <math.h>
#include <rsf.h>

#include "bilateral.h"

float bilateral(float *trace    /* input */,
		int ns          /* spray radius */,
		float ax        /* range radius */,
		float bx        /* domain radius */,
		int n           /* sample number */,
		int niter       /* number of iterations */,
		bool gauss      /* flag of Gaussian weight */)
/*< get a bilateral-weighted value >*/
{   
    int i1, i2, irep, is;
    float *output, t, norm;

    output = sf_floatalloc(2*ns+1);

    for (irep=0; irep < niter; irep++) {
/*	for (i1=0; i1 < n; i1++) {
	    output[i1] = trace[i1];
	}		
*/	/* loop over samples */
	for (i1=0; i1 < n; i1++) {	    
	    t = 0.;
	    norm = 0.;
	    /* loop over samples */
	    for (i2=-ns; i2<= ns; i2++) {
		if (i1+i2 >= 0 && i1+i2 < n) {
		    output[ns+i2]=trace[i1+i2];
		} else {
		    output[ns+i2]=0.;
		}
	    }	    
	    /* accumulate shifts */
	    for (is=-ns; is <= ns; is++) {
		if (gauss) {
		    t += output[ns+is]*expf(-0.5*is*is/(ax*ax+FLT_EPSILON))
			*expf(-0.5*(output[ns+is]-output[ns])
			      *(output[ns+is]-output[ns])/(bx*bx+FLT_EPSILON));
		    norm += expf(-0.5*is*is/(ax*ax+FLT_EPSILON))*
			expf(-0.5*(output[ns+is]-output[ns])
			     *(output[ns+is]-output[ns])/(bx*bx+FLT_EPSILON));
		} else {
		    t += output[ns+is]*(1.-fabsf(1.*is)/(ns+FLT_EPSILON))
			*expf(-bx*(output[ns+is]-output[ns])
			      *(output[ns+is]-output[ns]));
		    norm += (1.-fabsf(1.*is)/(ns+FLT_EPSILON))
			*expf(-bx*(output[ns+is]-output[ns])
			      *(output[ns+is]-output[ns]));
		}	
	    }
	    /* Normalize */
	    trace[i1] = t / (norm+FLT_EPSILON);
	}
    }
    return trace[(n-1)/2];
    
}

/* 	$Id$	 */


