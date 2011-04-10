/* Linear adaptive PEF operators */
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
#include "apefs.h"

static float* sfilt;
static float* nfilt;
static float* tmp;
static int sf1, sf2, nf1, nf2, n1, n2;

void apefs_init (float* ss, 
		 float* nn, 
		 int sff1, int sff2, 
		 int nd1, int nd2,
		 int nff1, int nff2) 
/*< initialize with a pointer to a matrix >*/
{
    sfilt = ss;
    nfilt = nn;
    sf1   = sff1;
    sf2   = sff2;
    nf1   = nff1;
    nf2   = nff2;
    n1    = nd1;
    n2    = nd2;
    tmp = sf_floatalloc(n1*n2);
}

void spefs_lop (bool adj, bool add, 
		int nx, int ny, float* smm, float* sdd) 
/*< linear signal PEF operator >*/
{
    int i, j, k, l;
    float**** filt, **model, **data;

    filt = sf_floatalloc4(sf1,sf2,n1,n2);
    data = sf_floatalloc2(n1,n2);
    model = sf_floatalloc2(n1,n2);

    sf_adjnull (adj,add,nx,ny,smm,sdd);

    for (l=0; l < n2; l++) {
	for (k=0; k < n1; k++) {
	    for (j=0; j < sf2; j++) {
		for (i=0; i < sf1; i++) {
		    filt[l][k][j][i] = 
			sfilt[l*n1*sf2*sf1+
			      k*sf2*sf1+
			      j*sf1+
			      i];
		}
	    }
	    model[l][k] = smm[l*n1+k];
	    data[l][k] = sdd[l*n1+k];
	}
    }
    for (l=0; l < n2; l++) {
	for (k=0; k < n1; k++) {
	    for (j=0; j < sf2; j++) {
		for (i=-sf1/2; i < (sf1+1)/2; i++) {
		    /* zero value boundary conditions */
		    if (l+j < 0 || l+j >= n2 || k+i < 0 || k+i >= n1) {
			continue; 
		    }
		    if (adj) {
			model[l+j][k+i] += filt[l][k][j][i+sf1/2]*data[l][k];
		    } else {
			data[l][k] += filt[l][k][j][i+sf1/2]*model[l+j][k+i];
		    }
		}
	    }
	}
    }
    for (l=0; l < n2; l++) {
	for (k=0; k < n1; k++) {
	    smm[l*n1+k] = model[l][k];
	    sdd[l*n1+k] = data[l][k];
	}
    }
    free (***filt); free(**filt); free(*filt); free(filt);
    free (*model); free (model);
    free (*data); free (data);
}

void npefs_lop (bool adj, bool add, 
		int nx, int ny, float* nmm, float* ndd) 
/*< linear cascading noise PEF operator >*/
{
    sf_chain(npefs_oper,npefs_oper,adj,add,nx,ny,nx,nmm,ndd,tmp);
}

void npefs_oper (bool adj, bool add, 
		 int nx, int ny, float* mm, float* dd) 
/*< linear noise PEF operator >*/
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
			nfilt[l*n1*nf2*nf1+
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

void apefs_close () 
/*< free filter memory >*/
{
    free (sfilt);
    free (nfilt);
    free (tmp);
}

/* 	$Id$	 */
