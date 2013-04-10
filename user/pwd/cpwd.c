/* Matrix definition for plane-wave destruction in the case of constant dips */
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
/*^*/

#include "cpwd.h"
#include "apfilt.h"

#ifndef _cpwd_h

typedef struct Cpwd *cpwd; /* abstract data type */
/*^*/

#endif

struct Cpwd {
    int n, na;
    float *a;
};

cpwd cpwd_init(int n1 /* trace length */, 
	       int nw /* filter order */)
/*< initialize >*/
{
    cpwd w;

    w = (cpwd) sf_alloc(1,sizeof(*w));
    w->n = n1;
    w->na = 2*nw+1;
 
    w->a = sf_floatalloc (w->na);
    
    apfilt_init(1);

    return w;
}

void cpwd_close (cpwd w)
/*< free allocated storage >*/
{
    apfilt_close();
    free (w->a);
    free (w);
}

float cpwd_define (bool adj    /* adjoint flag */, 
		   cpwd w       /* cpwd object */, 
		   float pp     /* slope */, 
		   float* offd /* defined off-diagonal */)
/*< fill the matrix (returns the diagonal) >*/
{
    float b[3], diag;
    
    passfilter (pp, b);
    
    if (adj) {
	w->a[0] = b[2];
	w->a[1] = b[1];	
	w->a[2] = b[0];
    } else {
	w->a[0] = b[0];
	w->a[1] = b[1];	
	w->a[2] = b[2];
    }
  
    diag = 
	w->a[0]*w->a[0] + 
	w->a[1]*w->a[1] + 
	w->a[2]*w->a[2];
    offd[0] += 
	w->a[0]*w->a[1] + 
	w->a[1]*w->a[2];
    offd[1] +=
	w->a[0]*w->a[2];
    
    return diag;
}

void cpwd_set (bool adj   /* adjoint flag */,
	       cpwd w      /* cpwd object */, 
	       float* inp /* input */, 
	       float* out /* output */, 
	       float* tmp /* temporary storage */)
/*< matrix multiplication >*/
{
    int i, n;
    float *a;

    n = w->n;
    a = w->a;

    if (adj) {
	for (i=0; i < n; i++) {
	    tmp[i]=0.;
	}
	for (i=1; i < n-1; i++) {
	    tmp[i-1] += a[0]*out[i];
	    tmp[i  ] += a[1]*out[i];
	    tmp[i+1] += a[2]*out[i];
	}
	/* zero value */
	tmp[0  ] += a[1]*out[0];
	tmp[1  ] += a[2]*out[0];
	tmp[n-2] += a[0]*out[n-1];
	tmp[n-1] += a[1]*out[n-1];
	
	for (i=0; i < n; i++) {
	    inp[i]=0.;
	}
	for (i=1; i < n-1; i++) {
	    inp[i-1] += a[0]*tmp[i];
	    inp[i  ] += a[1]*tmp[i];
	    inp[i+1] += a[2]*tmp[i];
	}
	/* zero value */
	inp[0  ] += a[1]*tmp[0];
	inp[1  ] += a[2]*tmp[0];
	inp[n-2] += a[0]*tmp[n-1];
	inp[n-1] += a[1]*tmp[n-1];
    } else {
	for (i=1; i < n-1; i++) {
	    tmp[i] = 
		a[0]*inp[i-1] +
		a[1]*inp[i  ] +
		a[2]*inp[i+1];
	}
	/* zero value */
	tmp[0] = 
	    a[1]*inp[0] +
	    a[2]*inp[1];
	tmp[n-1] = 
	    a[0]*inp[n-2] +
	    a[1]*inp[n-1];

	for (i=1; i < n-1; i++) {
	    out[i] = 
		a[0]*tmp[i-1] +
		a[1]*tmp[i  ] +
		a[2]*tmp[i+1];
	}
	/* zero value */
	out[0] = 
	    a[1]*tmp[0] +
	    a[2]*tmp[1];
	out[n-1] = 
	    a[0]*tmp[n-2] +
	    a[1]*tmp[n-1];
    }
}
