/* Fourier-domain zero-offset velocity continuation */
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

#include "fint1.h"

static int nt, nx, n2, n3, nw;
static float dt, dx, t0, *strace, o2, d2, dw;
static sf_complex *ctrace;
static fint1 str, istr;
static kiss_fftr_cfg forw, invs;


void velcon_init(int nt1   /* time samples */,
		 int nx1   /* space samples */,
		 float dt1 /* time sampling */,
		 float dx1 /* space sampling */,
		 float t1  /* time origin */,
		 int pad1, 
		 int pad2  /* padding in time */,
		 int next /* trace extension */)
/*< initialize >*/
{
    nt = nt1;
    nx = nx1;
    dt = dt1;
    dx = dx1;
    t0 = t1;

    n2 = pad1;
    n3 = pad2;

    /* t squared stretch */
    o2 = t0*t0;
    d2 = t0+(nt-1)*dt;
    d2 = (d2*d2 - o2)/(n2-1);

    str = fint1_init(next,nt,0);
    istr = fint1_init(next,n2,0);

    /* FFT in time */
    forw = kiss_fftr_alloc(n3,0,NULL,NULL);
    invs = kiss_fftr_alloc(n3,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("KISS FFT allocation error");

    nw = n3/2+1;
    dw = 2*SF_PI/(n3*d2);

    /* cosine FT in space */
    sf_cosft_init(nx);
    dx = SF_PI/(kiss_fft_next_fast_size(nx-1)*dx1);

    strace = sf_floatalloc(n3);
    ctrace = sf_complexalloc(nw);
}

void velcon_close(void)
/*< free allocated storage >*/
{
    free(ctrace);
    free(strace);
    sf_cosft_close();
    free(invs);
    free(forw);
    fint1_close(istr);
    fint1_close(str);
}

void velcon(float ** data /* [nx][nt] */,
	    float v0      /* initial velocity */,
	    float v1      /* new velocity */)
/*< apply velocity continuation >*/
{
    int it, ix, i2, iw;
    float *trace, k, w, t;
    sf_complex shift;

    for (it=0; it < nt; it++) {
	sf_cosft_frw (data[0], it, nt);
    }

    for (ix=0; ix < nx; ix++) {
	/* loop over wavenumbers */
	k = ix*dx;
	k *= k*(v0*v0 - v1*v1)/16;

	trace = data[ix];

	/* stretch t -> t^2 */

	for (it=0; it < nt; it++) {
	    trace[it] /= nt;
	}

	fint1_set(str,trace);

	for (i2=0; i2 < n2; i2++) {
	    t = o2+i2*d2;
	    t = sqrtf(t);
	    t = (t-t0)/dt;
	    it = t;
	    if (it >= 0 && it < nt) {
		strace[i2] = fint1_apply(str,it,t-it,false);
	    } else {
		strace[i2] = 0.;
	    }
	}

	for (i2=n2; i2 < n3; i2++) {
	    strace[i2] = 0.;
	}

	/* FFT */

	kiss_fftr(forw,strace, (kiss_fft_cpx *) ctrace);

	/* velocity continuation */

	ctrace[0]=sf_cmplx(0.,0.); /* dc */

	for (iw=1; iw < nw; iw++) {
	    w = iw*dw;
	    w = k/w;
	    shift = sf_cmplx(cosf(w)/n3,sinf(w)/n3);

#ifdef SF_HAS_COMPLEX_H
	    ctrace[iw] *= shift;
#else
	    ctrace[iw] = sf_cmul(ctrace[iw],shift);
#endif
	}

	/* Inverse FFT */

	kiss_fftri(invs,(const kiss_fft_cpx *) ctrace, strace);

	/* inverse stretch t^2->t */

	fint1_set(istr,strace);

	for (it=0; it < nt; it++) {
	    t = t0+it*dt;
	    t = t*t;
	    t = (t-o2)/d2;
	    i2 = t;
	    if (i2 >= 0 && i2 < n2) {
		trace[it] = fint1_apply(istr,i2,t-i2,false);
	    } else {
		trace[it] = 0.;
	    }
	}	
    } /* ix */

    /* inverse cosine FT in space */
    for (it=0; it < nt; it++) {
	sf_cosft_inv (data[0], it, nt);
    }
}
