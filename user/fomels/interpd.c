/* Trace interpolation with plane-wave destruction */
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

#include "interpd.h"
#include "banded.h"
#include "bandpass.h"
#include "burg.h"
#include "pwd.h"

static int n1, nb;
static bands slv;
static float **offd, eps;
static pwd w1, w2;

void interp_init (int n   /* trace length */, 
		  float e /* regularization */)
/*< initialize >*/
{
    const int nw=1;

    n1 = n;
    nb = 2*nw;

    eps = e;
  
    slv = banded_init (n1, nb);
    offd = sf_floatalloc2 (n1,nb);

    bandpass_init();
    w1 = pwd_init (n1, nw);
    w2 = pwd_init (n1, nw);
}

void interp_close (void)
/*< free allocated storage >*/
{
    banded_close (slv);
    free (*offd);
    free (offd);
    pwd_close (w1);
    pwd_close (w2);
    bandpass_close ();
}

void interp2(int n2      /* number of traces */, 
	     float** in  /* input [n2][n1] */, 
	     float** out /* output [2*n2-1][n1] */, 
	     float** pp  /* slope [n2][n1] */)
/*< interpolate >*/
{
    int i1, i2;
    float *diag, a;

    for (i2=0; i2 < n2-1; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    out[2*i2][i1] = in[i2][i1];
	}
	diag = out[2*i2+1];

	a = pef_burg2(n1,2,in+i2);
	pef_define (n1, a, eps, diag, offd);

	pwd_define (true,  w1, pp[i2  ], diag, offd);
	pwd_define (false, w2, pp[i2+1], diag, offd);

	banded_define (slv, diag, offd);

	pwd_set (false, w1, in[i2],   offd[0], diag);
	pwd_set (false, w2, in[i2+1], offd[1], diag);
	for (i1=0; i1 < n1; i1++) {
	    diag[i1] = offd[0][i1] + offd[1][i1];
	}

	banded_solve (slv, diag);
	bandpass (n1, diag); 
    }
    for (i1=0; i1 < n1; i1++) {
	out[2*n2-2][i1] = in[n2-1][i1];
    }
}

