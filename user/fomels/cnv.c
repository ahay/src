/* Matrix definition for convection */
/*
  Copyright (C) 2025 University of Texas at Austin
  
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

#include "cnv.h"

#ifndef _cnv_h

typedef struct Cnv *cnv; /* abstract data type */
/*^*/

#endif

struct Cnv {
    int n, na;
    float **a, *b;
};

cnv cnv_init(int n1 /* trace length */, 
	     int nw /* filter order */)
/*< initialize >*/
{
    cnv w;

    w = (cnv) sf_alloc(1,sizeof(*w));
    w->n = n1;
    w->na = 2*nw+1;
 
    w->a = sf_floatalloc2 (n1,w->na);
    w->b = sf_floatalloc (w->na);

    return w;
}

void cnv_close (cnv w)
/*< free allocated storage >*/
{
    free (w->a[0]);
    free (w->a);
    free (w->b);
    free (w);
}

void cnv_define (bool adj     /* adjoint flag */, 
		 cnv w        /* cnv object */, 
		 float** cc   /* convection */, 
		 float* diag  /* defined diagonal */, 
		 float** offd /* defined off-diagonal */)
/*< fill the matrix >*/
{
  int i, j, k, m, n, iw, nw;
    float am, aj;
    
    nw = (w->na-1)/2; 
    n = w->n;

    for (i=0; i < n; i++) {
      w->b[nw] = 1.0f;
      for (iw=1; iw <= nw; iw++) {
	w->b[nw+iw] = cc[2*iw-2][i];
	w->b[nw-iw] = cc[2*iw-1][i];
	w->b[nw] -= w->b[nw-iw] + w->b[nw+iw];
      }
	
      for (j=0; j < w->na; j++) {
	if (adj) {
	  w->a[j][i] = w->b[w->na-1-j];
	} else {
	  w->a[j][i] = w->b[j];
	}
      }
    }
    
    for (i=0; i < n; i++) {
	for (j=0; j < w->na; j++) {
	    k = i+j-nw;
	    if (k >= nw && k < n-nw) {
		aj = w->a[j][k];
		diag[i] += aj*aj;
	    }
	} 
	for (m=0; m < 2*nw; m++) {
	    for (j=m+1; j < w->na; j++) {
		k = i+j-nw;
		if (k >= nw && k < n-nw) {
		    aj = w->a[j][k];
		    am = w->a[j-m-1][k];
		    offd[m][i] += am*aj;
		}
	    }
	}
    }
}

void cnv_set (bool adj   /* adjoint flag */,
	      cnv w      /* cnv object */, 
	      float* inp /* input */, 
	      float* out /* output */, 
	      float* tmp /* temporary storage */)
/*< matrix multiplication >*/
{
    int i, j, k, n, nw;

    nw = (w->na-1)/2;
    n = w->n;

    if (adj) {
	for (i=0; i < n; i++) {
	    tmp[i]=0.;
	}
	for (i=0; i < n; i++) {
	    for (j=0; j < w->na; j++) {
		k = i+j-nw;
		if (k >= nw && k < n-nw) 
		    tmp[k] += w->a[j][k]*out[i];
	    }
	}
	for (i=0; i < n; i++) {
	    inp[i]=0.;
	}
	for (i=nw; i < n-nw; i++) {
	    for (j=0; j < w->na; j++) {
		k = i+j-nw;
		inp[k] += w->a[j][i]*tmp[i];
	    }
	}
    } else {
	for (i=0; i < n; i++) {
	    tmp[i] = 0.;
	}
	for (i=nw; i < n-nw; i++) {
	    for (j=0; j < w->na; j++) {
		k = i+j-nw;
		tmp[i] += w->a[j][i]*inp[k];
	    }
	}
	for (i=0; i < n; i++) {
	    out[i] = 0.;
	    for (j=0; j < w->na; j++) {
		k = i+j-nw;
		if (k >= nw && k < n-nw) 
		    out[i] += w->a[j][k]*tmp[k];
	    }
	}
    }
}
