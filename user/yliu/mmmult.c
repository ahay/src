/* Multi-dimension matrix multiplication operator */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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
#include "mmmult.h"

static float* filter;
static int nf1, nf2, n1, n2;

void mmmult_init (float* bb, int nff1, int nff2, int nff3, int nff4) 
/*< initialize with a pointer to a matrix >*/
{
    filter = bb;
    nf1    = nff1;
    nf2    = nff2;
    n1     = nff3;
    n2     = nff4;
}

void mmmult_lop (bool adj, bool add, 
		  int nx, int ny, float* mm, float* dd) 
/*< linear operator >*/
{
    int i, j, k, l;
    float**** filt, **model, **data;

    filt = sf_floatalloc4(nf1,nf2,n1,n2);
    data = sf_floatalloc2(n1,n2);
    model = sf_floatalloc2(n1,n2);

    sf_adjnull (adj,add,nx,ny,mm,dd);

    for (l=0; l < n2; l++) {
	for (k=0; k < n1; k++) {
	    for (j=0; j < nf2; j++) {
		for (i=0; i < nf1; i++) {
		    filt[l][k][j][i] = 
			filter[l*n1*nf2*nf1+
			       k*nf2*nf1+
			       j*nf1+
			       i];
		}
	    }
	    model[l][k] = mm[l*n1+k];
	    data[l][k] = dd[l*n1+k];
	}
    }
    for (l=0; l < n2; l++) {
	for (k=0; k < n1; k++) {
	    for (j=0; j < nf2; j++) {
		for (i=-nf1/2; i < (nf1+1)/2; i++) {
		    /* zero value boundary conditions */
		    if (l+j < 0 || l+j >= n2 || k+i < 0 || k+i >= n1) {
			continue; 
		    }
		    if (adj) {
			model[l+j][k+i] += filt[l][k][j][i+nf1/2]*data[l][k];
		    } else {
			data[l][k] += filt[l][k][j][i+nf1/2]*model[l+j][k+i];
		    }
		}
	    }
	}
    }
    for (l=0; l < n2; l++) {
	for (k=0; k < n1; k++) {
	    mm[l*n1+k] = model[l][k];
	    dd[l*n1+k] = data[l][k];
	}
    }
    free (***filt); free(**filt); free(*filt); free(filt);
    free (*model); free (model);
    free (*data); free (data);
}

void mmmult_close () 
/*< free filter memory >*/
{
    free (filter);
}

/* 	$Id: mmmult.c 7107 2011-04-10 02:04:14Z ivlad $	 */
