/* Trace prediction with plane-wave destruction */
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

#include "predict.h"
#include "banded.h"
#include "pwd.h"

static int n1, n2, nb;
static bands slv;
static float *diag, **offd, eps, **dip, *tt;
static pwd w;

void predict_init (int nx, int ny /* data size */, 
		   float e        /* regularization parameter */)
/*< initialize >*/
{
    const int nw=1;

    n1 = nx;
    n2 = ny;
    nb = 2*nw;

    eps = e;

    slv = banded_init (n1, nb);
    diag = sf_floatalloc (n1);
    offd = sf_floatalloc2 (n1,nb);

    w = pwd_init (n1, nw);
    tt = NULL;
}

void predict_close (void)
/*< free allocated storage >*/
{
    banded_close (slv);
    free (diag);
    free (*offd);
    free (offd);
    pwd_close (w);
    if (NULL != tt) free(tt);
}

void predict_step(bool adj     /* adjoint flag */,
		  bool forw    /* forward or backward */, 
		  float* trace /* trace */, 
		  float* pp    /* slope */)
/*< prediction step >*/
{
    int i1;

    for (i1=0; i1 < n1; i1++) {
	diag[i1] = 6.*eps;
	offd[0][i1] = -4.*eps;
	offd[1][i1] = eps;
    }
    diag[0] = diag[n1-1] = 1.+eps;
    diag[1] = diag[n1-2] = 1.+5.*eps;
    offd[0][0] = offd[0][n1-2] = -2.*eps;

    pwd_define (forw, w, pp, diag, offd);
    banded_define (slv, diag, offd);

    if (adj) {
	banded_solve (slv, trace);
	pwd_set (true, w, trace, trace, diag);
    } else {
	pwd_set (false, w, trace, trace, diag);
	banded_solve (slv, trace);
    }
}

void predict_set(float **dip1 /* dip field [n2][n1] */)
/*< set the local slopes for applying the linear operator >*/
{
    dip=dip1;
    tt = sf_floatalloc(n1);
}

void predict_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy)
/*< linear operator >*/
{
    int i1, i2;

    if (nx != ny || nx != n1*n2) sf_error("%s: Wrong dimensions",__FILE__);

    sf_adjnull(adj,add,nx,ny,xx,yy);

    for (i1=0; i1 < n1; i1++) {
	tt[i1] = 0.;
    }

    if (adj) {
	for (i2=n2-1; i2 >= 0; i2--) {
	    predict_step(true,true,tt,dip[i2]);
	    for (i1=0; i1 < n1; i1++) {
		tt[i1] += yy[i1+i2*n1];
	    }
	    for (i1=0; i1 < n1; i1++) {
		xx[i1+i2*n1] += tt[i1];
	    }
	}
    } else {
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		tt[i1] += xx[i1+i2*n1];
	    }
	    for (i1=0; i1 < n1; i1++) {
		yy[i1+i2*n1] += tt[i1];
	    }
	    predict_step(false,true,tt,dip[i2]);
	}
    }
}

void predict_flat(int i0     /* reference trace number */, 
		  float** d  /* input */, 
		  float** m  /* output */, 
		  float** pp /* slope */)
/*< predictive flattening >*/
{
    int i1, i2, k2;
    float *trace;

    /* prediction from the left */
    for (i2=0; i2 <= i0; i2++) {
        for (i1=0; i1 < n1; i1++) {
            m[i2][i1] = d[i2][i1];
        }

        if (i2 == i0) break;

	for (i1=0; i1 < n1; i1++) {
	    diag[i1] = 6.*eps;
	    offd[0][i1] = -4.*eps;
	    offd[1][i1] = eps;
	}
	diag[0] = diag[n1-1] = 1.+eps;
	diag[1] = diag[n1-2] = 1.+5.*eps;
	offd[0][0] = offd[0][n1-2] = -2.*eps;

        pwd_define (true, w, pp[i2], diag, offd);
        banded_define (slv, diag, offd);

        for (k2=0; k2 <= i2; k2++) {
            trace = m[k2];

            pwd_set (false, w, trace, trace, diag);
            banded_solve (slv, trace);
        }
    }
    
    /* prediction from the right */
    for (i2=n2-1; i2 > i0; i2--) {
        for (i1=0; i1 < n1; i1++) {
            m[i2][i1] = d[i2][i1];
        }

	for (i1=0; i1 < n1; i1++) {
	    diag[i1] = 6.*eps;
	    offd[0][i1] = -4.*eps;
	    offd[1][i1] = eps;
	}
	diag[0] = diag[n1-1] = 1.+eps;
	diag[1] = diag[n1-2] = 1.+5.*eps;
	offd[0][0] = offd[0][n1-2] = -2.*eps;

        pwd_define (false, w, pp[i2-1], diag, offd);
        banded_define (slv, diag, offd);

        for (k2=n2-1; k2 >= i2; k2--) {
            trace = m[k2];

            pwd_set (false, w, trace, trace, diag);
            banded_solve (slv, trace);
        }
    }
}
