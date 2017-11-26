/* N-D triangle smoothing derivative */
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

#include "dtrianglen.h"
#include "file.h"
#include "triangle.h"
#include "decart.h"
#include "alloc.h"
#include "deriv.h"

static int *n, s[SF_MAX_DIM], *rect, nd, dim;
static sf_triangle *tr;

void sf_dtrianglen_init (int ndim  /* number of dimensions */, 
			 int *nbox /* triangle radius [ndim] */, 
			 int *ndat /* data dimensions [ndim] */)
/*< initialize >*/
{
    int i;

    dim = ndim;
    n = sf_intalloc(dim);

    tr = (sf_triangle*) sf_alloc(dim,sizeof(sf_triangle));
    rect = nbox;

    nd = 1;
    for (i=0; i < dim; i++) {
	tr[i] = (nbox[i] > 1)? sf_triangle_init (nbox[i],ndat[i],false): NULL;
	s[i] = nd;
	n[i] = ndat[i];
	nd *= ndat[i];
    }
}

void sf_dtrianglen (int ider   /* direction of the derivative */,
		    int nrep   /* how many times to repeat smoothing */,
		    int nderiv /* derivative filter accuracy */,
		    float* data   /* input/output */)
/*< linear operator (derivative with respect to radius) >*/
{
    float *t1=NULL, *t2=NULL;
    int i, irep, j, i0, i1, n1=0, s1=0;

    ider--;

    if (ider >= 0) {
	n1 = n[ider];
	s1 = s[ider];

	t1 = sf_floatalloc(n1);
	t2 = sf_floatalloc(n1);

	sf_deriv_init(n1,nderiv,0.);
    } 

    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) {
	    for (j=0; j < nd/n[i]; j++) {
		i0 = sf_first_index (i,j,dim,n,s);
		if (i==ider) {
		    for (i1=0; i1 < n1; i1++) {
			t1[i1] = data[i0+i1*s1];
		    }
		    for (irep=0; irep < nrep-1; irep++) {
			sf_smooth2 (tr[i], 0, 1, false, t1);
		    }
		    sf_dsmooth (tr[i],0,1,false,t1);
		    sf_deriv(t1,t2);
		}
		
		for (irep=0; irep < nrep; irep++) {
		    sf_smooth2 (tr[i], i0, s[i], false, data);
		}
		
		if (i==ider) {
		    for (i1=0; i1 < n1; i1++) {
			data[i0+i1*s1] = nrep*(t2[i1] - 2*data[i0+i1*s1]/rect[i]);
		    }
		}
	    }
	}
    }

    if (ider >= 0) {
	free(t1);
	free(t2);
	sf_deriv_close();
    }
}

void sf_dtrianglen_close(void)
/*< free allocated storage >*/
{
    int i;

    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) sf_triangle_close (tr[i]);
    }

    free(tr);
    free(n);
}

