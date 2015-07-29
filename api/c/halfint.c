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
#include <math.h>

#include "halfint.h"
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include "alloc.h"
#include "error.h"
#include "_defs.h"
#include "komplex.h"
#include "adjnull.h"

#include "_bool.h"
/*^*/

static int nw, n;
static float *tmp;
static kiss_fft_cpx *cx, *cf;
static kiss_fftr_cfg forw, invs;

void sf_halfint_init (bool inv  /* differentiation or integration */, 
		      int n1    /* trace length */, 
		      float rho /* regularization */)
/*< Initialize >*/
{
    int i;
    float om;
    kiss_fft_cpx cw, cz, cz2;

    n = 2*kiss_fft_next_fast_size((n1+1)/2);
    nw = n/2+1;

    cx = (kiss_fft_cpx*) sf_complexalloc(nw);
    cf = (kiss_fft_cpx*) sf_complexalloc(nw);

    forw = kiss_fftr_alloc(n,0,NULL,NULL);
    invs = kiss_fftr_alloc(n,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error");

    for (i=0; i < nw; i++) {
	om = -2.*SF_PI*i/n;
	cw.r = cosf(om);
	cw.i = sinf(om);

	cz.r = 1.-rho*cw.r;
	cz.i = -rho*cw.i;
	if (inv) {
	    cz = sf_csqrtf(cz);
	} else {
	    cz2.r = 0.5*(1.+rho*cw.r);
	    cz2.i = 0.5*rho*cw.i;
	    cz = sf_csqrtf(sf_cdiv(cz2,cz));
	}
	cf[i].r = cz.r/n;
	cf[i].i = cz.i/n;
    }

    tmp = sf_floatalloc(n);
}

void sf_halfint_lop(bool adj, bool add, int n1, int n2, float *xx, float *yy)
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

    sf_halfint (adj,tmp);

    for (i=0; i < n1; i++) {
	if (adj) {
	    xx[i] += tmp[i];
	} else {
	    yy[i] += tmp[i];
	}
    }
}

void sf_halfint (bool adj, float* x /* [n] */)
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

void sf_halfint_close(void)
/*< Free allocated storage >*/
{
    free (cx);
    free (cf);
    free (forw);
    free (invs);
    free (tmp);
}
