/* N-D non-stationary triangle smoothing as a linear operator */
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

#include "cntrianglen.h"

#include <rsf.h>
/*^*/

#include "cntriangle.h"

static int *n, s[SF_MAX_DIM], nd, dim, **tsft, nrep;
static cntriangle *tr;
static sf_complex *tmp;
static float **tlen;

void cntrianglen_init (int ndim  /* number of dimensions */, 
		      int *nbox /* triangle radius [ndim] */, 
		      int *ndat /* data dimensions [ndim] */,
		      float **len /* triangle lengths [ndim][nd] */,
                      int **sft /* triangle shifts [ndim][nd] */,
		      int repeat /* repeated smoothing */)
/*< initialize >*/
{
    int i;

    n = ndat;
    dim = ndim;

    tr = (cntriangle*) sf_alloc(dim,sizeof(cntriangle));
    
    nd = 1;
    for (i=0; i < dim; i++) {
	tr[i] = (nbox[i] > 1)? cntriangle_init (nbox[i],ndat[i]): NULL;
	s[i] = nd;
	nd *= ndat[i];
    }
    tlen = len; 
    tsft = sft;

    tmp = sf_complexalloc(nd);
    nrep = repeat;
}

void cntrianglen_lop (bool adj, bool add, int nx, int ny, sf_complex* x, sf_complex* y)
/*< linear operator >*/
{
    int i, j, i0, irep;

    if (nx != ny || nx != nd) 
	sf_error("%s: Wrong data dimensions: nx=%d, ny=%d, nd=%d",
		 __FILE__,nx,ny,nd);

    sf_cadjnull (adj,add,nx,ny,x,y);
  
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
	    for (j=0; j < nd/n[i]; j++) {
		i0 = sf_first_index (i,j,dim,n,s);

		for (irep=0; irep < nrep; irep++) {
		    cnsmooth (tr[i], i0, s[i], false, tlen[i], tsft[i], tmp);
		}
	    }
	}
    }
	
    if (adj) {
	for (i=0; i < nd; i++) {
#ifdef SF_HAS_COMPLEX_H
	    x[i] += tmp[i];
#else
	    x[i] = sf_cadd(x[i],tmp[i]);
#endif
	}
    } else {
	for (i=0; i < nd; i++) {
#ifdef SF_HAS_COMPLEX_H
	    y[i] += tmp[i];
#else
	    y[i] = sf_cadd(y[i],tmp[i]);
#endif
	}
    }     
}

void cntrianglen_close(void)
/*< free allocated storage >*/
{
    int i;

    free (tmp);

    for (i=0; i < dim; i++) {
	if (NULL != tr[i]) cntriangle_close (tr[i]);
    }

    free(tr);
}

