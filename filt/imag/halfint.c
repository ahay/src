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

static int nn, nw;
static float complex *cx, *cf;
static kiss_fftr_cfg forw, invs;

void halfint_init (bool adj  /* causal or anticausal */, 
		   bool inv  /* differentiation or integration */, 
		   int n     /* trace length */, 
		   float rho /* regularization */)
/*< Initialize >*/
{
    int i;
    float om;
    float complex cz;

    nn = n;
    nw = nn/2+1;

    cx = sf_complexalloc(nw);
    cf = sf_complexalloc(nw);

    forw = kiss_fftr_alloc(nn,0,NULL,NULL);
    invs = kiss_fftr_alloc(nn,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("%s: KISS FFT allocation error");

    for (i=0; i < nw; i++) {
	om = 2.*SF_PI*i/nn;
        if (!adj) om = - om;

	cz = cexpf(I*om);
	if (inv) {
	    cf[i] = csqrtf(1.-rho*cz);
	} else {
	    cf[i] = csqrtf(0.5*(1.+rho*cz)/(1.-rho*cz));
	}
    }
}

void halfint (float* x /* [n] */)
/* < Integrate in place >*/
{
    int i;

    kiss_fftr(forw,x, (kiss_fft_cpx *) cx);
    for (i=0; i < nw; i++) {
	cx[i] *= cf[i]/nn;
    }
    kiss_fftri(invs,(const kiss_fft_cpx *) cx, x);
}

void halfint_close(void)
/*< Free allocated storage >*/
{
    free (cx);
    free (cf);
    free (forw);
    free (invs);
}

