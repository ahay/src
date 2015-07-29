/* N-D triangle smoothing as a linear operator */
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

#include "trianglen.h"
#include "triangle.h"

static int *n, s[SF_MAX_DIM], nd, *nb, dim;
static int *nn, **i0, **il;
static float **tr;
static float *tmp;

void trianglen_init (int ndim  /* number of dimensions */, 
		     int *nbox /* triangle radius [ndim] */, 
		     int *ndat /* data dimensions [ndim] */)
/*< initialize >*/
{
    int i;

    dim = ndim;
    nb  = nbox;
    n = sf_intalloc(dim);    

    tr = (float **) sf_alloc(dim,sizeof(float *));

    nd = 1;
    for (i=0; i < dim; i++) {
	tr[i] = (nbox[i] > 1)? sf_floatalloc (ndat[i]+2*nbox[i]): NULL;
	s[i] = nd;
	n[i] = ndat[i];
	nd *= ndat[i];
    }
    tmp = sf_floatalloc (nd);
}

void trianglen_topo (int *topo /* topography */,
		     int seg   /* maximum number of segments */)
/*< pre-process topography >*/
{
    int i, j, k, t, size;
    bool new;

    nn = sf_intalloc (dim);

    size = 1;
    for (i=0; i < dim; i++) {
	size = (nd/n[i] > size)? (nd/n[i]): size;
    }
    
    i0 = sf_intalloc2 (size*seg,dim);
    il = sf_intalloc2 (size*seg,dim);

    for (i=0; i < dim; i++) {
	nn[i] = 0;
	
	if (NULL != tr[i]) {
	    for (j=0; j < nd/n[i]; j++) {
		k = sf_first_index (i,j,dim,n,s);
		
		new = true;
		for (t=0; t < n[i]; t++) {
		    if (topo[k+t*s[i]]==1 && new) {
			i0[i][nn[i]] = k+t*s[i];
			il[i][nn[i]] = 1;
			nn[i]++;
			new = false;
			continue;
		    }
		    if (topo[k+t*s[i]]!=1) new = true;
		    if (!new) il[i][nn[i]-1]++;
		}
	    }
	}
    }
}

void trianglen_lop (bool adj, bool add, int nx, int ny, float* x, float* y)
/*< linear operator >*/
{
    int i, j;

    if (nx != ny || nx != nd) 
	sf_error("%s: Wrong data dimensions: nx=%d, ny=%d, nd=%d",
		 __FILE__,nx,ny,nd);

    sf_adjnull (adj,add,nx,ny,x,y);
  
    if (adj) {
	for (i=0; i < nd; i++) {
	    tmp[i] = y[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    tmp[i] = x[i];
	}
    }

  
    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) {
	    for (j=0; j < nn[i]; j++) {
		smooth (tr[i],i0[i][j],s[i],il[i][j],nb[i],false,false,tmp);
	    }
	}
    }
	
    if (adj) {
	for (i=0; i < nd; i++) {
	    x[i] += tmp[i];
	}
    } else {
	for (i=0; i < nd; i++) {
	    y[i] += tmp[i];
	}
    }    
}

void trianglen_close(void)
/*< free allocated storage >*/
{
    int i;

    free (tmp);

    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) free (tr[i]);
	free (i0[i]);
	free (il[i]);
    }

    free(tr);
    free(n);
    free(i0);
    free(il);
}

/* 	$Id: trianglen.c 7107 2011-04-10 02:04:14Z ivlad $	 */
