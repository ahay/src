/* Matrix definition for plane-wave destruction */
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

#include "pwd.h"
#include "apfilt.h"

#ifndef _pwd_h

typedef struct Pwd *pwd; /* abstract data type */
/*^*/

#endif

struct Pwd {
    int n, na;
    float **a, *b;
};

pwd pwd_init(int n1 /* trace length */, 
	     int nw /* filter order */)
/*< initialize >*/
{
    pwd w;

    w = (pwd) sf_alloc(1,sizeof(*w));
    w->n = n1;
    w->na = 2*nw+1;
 
    w->a = sf_floatalloc2 (n1,w->na);
    w->b = sf_floatalloc (w->na);

    apfilt_init(nw);

    return w;
}

void pwd_close (pwd w)
/*< free allocated storage >*/
{
    free (w->a[0]);
    free (w->a);
    free (w->b);
    free (w);
    apfilt_close();
}

void pwd_define (bool adj        /* adjoint flag */, 
		 pwd w           /* pwd object */, 
		 const float* pp /* slope */, 
		 float* diag     /* defined diagonal */, 
		 float** offd    /* defined off-diagonal */)
/*< fill the matrix >*/
{
    int i, j, k, nw;
    
    nw = (w->na-1)/2;

    for (i=0; i < w->n; i++) {
	passfilter (pp[i], w->b);
	
	for (j=0; j < w->na; j++) {
	    if (adj) {
		w->a[j][i] = w->b[w->na-1-j];
	    } else {
		w->a[j][i] = w->b[j];
	    }
	}
    }
    
    for (i=nw; i < w->n-nw; i++) {
	for (j=0; j < w->na; j++) {
	    diag[i] += w->a[j][i-nw+j]*w->a[j][i-nw+j]; 
	}
	for (k=0; k < 2*nw; k++) {
	    for (j=0; j < 2*nw-k; j++) {
		offd[k][i-k] += w->a[j][i+j]*w->a[j+k+1][i+j];
	    }
	}
/*
	offd[0][i] += 
	    w->a[0][i  ]*w->a[1][i  ] + 
	    w->a[1][i+1]*w->a[2][i+1];
	offd[1][i-1] +=
	    w->a[0][i]*w->a[2][i];
*/
    }
    /* zero-slope */
    offd[0][0] += w->a[0][0]*(w->a[1][0]+w->a[2][0]) + w->a[1][1]*w->a[2][1];
    offd[0][w->n-1] += w->a[0][w->n-1]*w->a[2][w->n-1];
    diag[0] += 
	(w->a[1][0]+w->a[2][0])*(w->a[1][0]+w->a[2][0]) + 
	w->a[2][1]*w->a[2][1];
    diag[w->n-1] += 
	w->a[0][w->n-2]*w->a[0][w->n-2] + 
	(w->a[0][w->n-1] + w->a[1][w->n-1])*(w->a[0][w->n-1] + w->a[1][w->n-1]);
    /* zero value 
    offd[0][0] += 
	w->a[0][0]*w->a[1][0] + 
	w->a[1][1]*w->a[2][1];
    diag[0] += 
	w->a[1][0]*w->a[1][0] + 
	w->a[2][1]*w->a[2][1];
    diag[w->n-1] += 
	w->a[0][w->n-2]*w->a[0][w->n-2] + 
	w->a[1][w->n-1]*w->a[1][w->n-1];
    */
}

void pwd_set (bool adj   /* adjoint flag */,
	      pwd w      /* pwd object */, 
	      float* inp /* input */, 
	      float* out /* output */, 
	      float* tmp /* temporary storage */)
/*< matrix multiplication >*/
{
    int i, j, k, nw;

    nw = (w->na-1)/2;

    if (adj) {
	for (i=0; i < w->n; i++) {
	    tmp[i]=0.;
	}
	for (i=0; i < w->n; i++) {
	    for (j=0; j < w->na; j++) {
		k = i+j-nw;
		if (k < 0) {
		    k = -k-1;
		} else if (k >= w->n) {
		    k = 2*w->n-k-1;
		}
		tmp[k] += w->a[j][k]*out[i];
	    }
	}
	for (i=0; i < w->n; i++) {
	    inp[i]=0.;
	}
	for (i=0; i < w->n; i++) {
	    for (j=0; j < w->na; j++) {
		k = i+j-nw;
		if (k < 0) {
		    k = -k-1;
		} else if (k >= w->n) {
		    k = 2*w->n-k-1;
		}
		inp[k] += w->a[j][i]*tmp[i];
	    }
	}
    } else {
	for (i=0; i < w->n; i++) {
	    tmp[i] = 0.;
	    for (j=0; j < w->na; j++) {
		k = i+j-nw;
		if (k < 0) {
		    k = -k-1;
		} else if (k >= w->n) {
		    k = 2*w->n-k-1;
		}
		tmp[i] += w->a[j][i]*inp[k];
	    }
	}
	for (i=0; i < w->n; i++) {
	    out[i] = 0.;
	    for (j=0; j < w->na; j++) {
		k = i+j-nw;
		if (k < 0) {
		    k = -k-1;
		} else if (k >= w->n) {
		    k = 2*w->n-k-1;
		}
		out[i] += w->a[j][k]*tmp[k];
	    }
	}
    }
}
