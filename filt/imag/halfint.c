/* Half-order integration and differentiation. */
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

#include "halfint.h"

#include <math.h>

#include <rsf.h>
/*^*/

static int nw, n;
static float *tmp;
static kiss_fft_cpx *cx, *cf;
static kiss_fftr_cfg forw, invs;

void halfint_init (bool inv  /* differentiation or integration */, 
		   int n1    /* trace length */, 
		   int m     /* order */,
		   float rho /* regularization */)
/*< Initialize >*/
{
    int i, j;
    float om, z;
    kiss_fft_cpx cw, cz, cz2;

    n = n1;
    nw = n/2+1;

    cx = (kiss_fft_cpx*) sf_complexalloc(nw);
    cf = (kiss_fft_cpx*) sf_complexalloc(nw);

    forw = kiss_fftr_alloc(n,0,NULL,NULL);
    invs = kiss_fftr_alloc(n,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error");

    for (i=0; i < nw; i++) {
	om = 2.*SF_PI*i/n;
	cw.r = cosf(om);
	cw.i = sinf(om);

	if (inv) {
	    z = 1.;
	    for (j=m; j >= 1; j--) {
		z = 1.+z*j/(2*j+1)*(1-cw.r);
	    }
	    cz.r=0;
	    cz.i=cw.i*z;
	    cz = sf_csqrtf(cz);
	} else {
	    cz.r = 1.-rho*cw.r;
	    cz.i = rho*cw.i;
	    cz2.r = 0.5*(1.+rho*cw.r);
	    cz2.i = -0.5*rho*cw.i;
	    cz = sf_csqrtf(sf_cdiv(cz2,cz));
	}
	cf[i].r = cz.r/n;
	cf[i].i = cz.i/n;
    }

    tmp = sf_floatalloc(n);
}

void halfint_lop(bool adj, bool add, int n1, int n2, float *xx, float *yy)
/*< linear operator >*/
{
    int i;

    sf_adjnull (adj,add,n1,n2,xx,yy);
    
    for (i=0; i < n1; i++) {
	tmp[i] = adj? yy[i]: xx[i];
    }
    for (i=n1; i < n; i++) {
	tmp[i] = 0.;
    }

    halfint (adj,tmp);

    for (i=0; i < n1; i++) {
	if (adj) {
	    xx[i] += tmp[i];
	} else {
	    yy[i] += tmp[i];
	}
    }
}

void halfint (bool adj, float* x /* [n] */)
/*< Integrate in place >*/
{
    int i;

    kiss_fftr(forw, x, cx);
    for (i=0; i < nw; i++) {
	if (adj) {
	    cx[i] = sf_cmul(cx[i],sf_conjf(cf[i]));
	} else {
	    cx[i] = sf_cmul(cx[i],cf[i]);
	}
    }
    kiss_fftri(invs, cx, x);
}

void halfint_close(void)
/*< Free allocated storage >*/
{
    free (cx);
    free (cf);
    free (forw);
    free (invs);
    free (tmp);
}
