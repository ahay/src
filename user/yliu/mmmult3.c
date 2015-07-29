/* 3-D matrix multiplication operator */
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
#include "mmmult3.h"

static float* filter;
static int nf1, nf2, nf3, n1, n2, n3;

void mmmult3_init (float* bb, 
		  int nff1, int nff2, int nff3, 
		  int nff4, int nff5, int nff6) 
/*< initialize with a pointer to a matrix >*/
{
    filter = bb;
    nf1    = nff1;
    nf2    = nff2;
    nf3    = nff3;
    n1     = nff4;
    n2     = nff5;
    n3     = nff6;
}

void mmmult3_lop (bool adj, bool add, 
		  int nx, int ny, float* mm, float* dd) 
/*< linear operator >*/
{
    int i, j, k, l, m, n;
    float****** filt, ***model, ***data;

    filt = sf_floatalloc6(nf1,nf2,nf3,n1,n2,n3);
    data = sf_floatalloc3(n1,n2,n3);
    model = sf_floatalloc3(n1,n2,n3);

    sf_adjnull (adj,add,nx,ny,mm,dd);

    for (n=0; n < n3; n++) {
	for (l=0; l < n2; l++) {
	    for (k=0; k < n1; k++) {
		for (m=0; m < nf3; m++) {
		    for (j=0; j < nf2; j++) {
			for (i=0; i < nf1; i++) {
			    filt[n][l][k][m][j][i] = 
				filter[n*n2*n1*nf3*nf2*nf1+
				       l*n1*nf3*nf2*nf1+
				       k*nf3*nf2*nf1+
				       m*nf2*nf1+
				       j*nf1+
				       i];
			}
		    }
		}
		model[n][l][k] = mm[n*n2*n1+l*n1+k];
		data[n][l][k] = dd[n*n2*n1+l*n1+k];
	    }
	}
    }
    for (n=0; n < n3; n++) {
	for (l=0; l < n2; l++) {
	    for (k=0; k < n1; k++) {
		for (m=0; m < nf3; m++) {
		    for (j=0; j < nf2; j++) {
			for (i=-nf1/2; i < (nf1+1)/2; i++) {
			    /* zero value boundary conditions */
			    if (l+j < 0 || l+j >= n2 || 
				k+i < 0 || k+i >= n1 ||
				n+m < 0 || n+m >= n3) {
				continue; 
			    }
			    if (adj) {
				model[n+m][l+j][k+i] += 
				    filt[n][l][k][m][j][i+nf1/2]*
				    data[n][l][k];
			    } else {
				data[n][l][k] += 
				    filt[n][l][k][m][j][i+nf1/2]*
				    model[n+m][l+j][k+i];
			    }
			}
		    }
		}
	    }
	}
    }
    for (n=0; n < n3; n++) {
	for (l=0; l < n2; l++) {
	    for (k=0; k < n1; k++) {
		mm[n*n2*n1+l*n1+k] = model[n][l][k];
		dd[n*n2*n1+l*n1+k] = data[n][l][k];
	    }
	}
    }
    free (*****filt); free (****filt); free(***filt); 
    free (**filt);    free (*filt);    free (filt);
    free (**model);   free (*model);   free (model);
    free (**data);    free (*data);    free (data);
}

void mmmult_close () 
/*< free filter memory >*/
{
    free (filter);
}

/* 	$Id: mmmult3.c 7107 2011-04-10 02:04:14Z ivlad $	 */
