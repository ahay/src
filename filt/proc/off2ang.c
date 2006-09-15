/* Offset to angle transformation for common-azimuth migration */
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

#include "off2ang.h"

static int nt, nt2, np, nx, nw;
static float *tt, **pp, dw, dp, p0, dx, x0;
static sf_complex **cd, **cm, *dd;
static kiss_fftr_cfg forw, invs;

void off2ang_init(int   nz /* depth samples */,
		  float dz /* depth sampling */,
		  int   na /* angle samples */,
                  float da /* angle sampling (in radians) */,
		  float a0 /* first angle (in radians) */,
		  int   nh /* offset sampling */,
		  float dh /* offset sampling */,
		  float h0 /* first offset */)
/*< initialize >*/
{
    nt=nz;
    np=na; dp=da; p0=a0;
    nx=nh; dx=dh; x0=h0;
    nt2 = 2*kiss_fft_next_fast_size(nt); /* padding */
    nw = nt2/2+1;
    dw = 2.0*SF_PI/(nt2*dz);
    
    tt = sf_floatalloc (nt2);
    pp = sf_floatalloc2 (np,nt);

    forw = kiss_fftr_alloc(nt2,0,NULL,NULL);
    invs = kiss_fftr_alloc(nt2,1,NULL,NULL);

    cm = sf_complexalloc2 (nw,np);
    cd = sf_complexalloc2 (nw,nx);
    dd = sf_complexalloc(nx);
}

void off2ang_close()
/*< free allocated storage >*/
{
    free(tt);
    free(*pp);
    free(pp);
    free(forw);
    free(invs);
    free(*cm);
    free(cm);
    free(*cd);
    free(cd);
    free(dd);
}

void off2ang(float **off /* input: local offset gather [nx][nz] */, 
	     float **ang /* output: angle gather [na][nz] */,
	     const float *dip  /* cross-line dip [nz] */)
/*< transform >*/
{
    int ix, ip, it, iw, iq;
    float w, p, dq, t;
    sf_complex m, c;

    /*** frequency-domain Radon transform ***/
    for (ix=0; ix < nx; ix++) { /* loop over offsets */	
	for (it=0; it < nt; it++) {
	    tt[it]=off[ix][it];
	}
	for (it=nt; it < nt2; it++) {
	    tt[it]=0.;
	}
		
	/* FFT to frequency */
	kiss_fftr(forw,tt, (kiss_fft_cpx *) cd[ix]);
    }

    for (iw=0; iw < nw; iw++) { /* loop over frequencies */
	w = iw*dw;

	for (ix=0; ix < nx; ix++) { /* loop over offsets */
            /* transpose + FFT scaling */
#ifdef SF_HAS_COMPLEX_H
	    dd[ix] = cd[ix][iw]/nt2;
#else
	    dd[ix] = sf_crmul(cd[ix][iw],1./nt2);
#endif
	}

	for (ip=0; ip < np; ip++) { /* loop over slopes */
	    p = tanf(p0+ip*dp);
	    t = w*p*dx;
	    c = sf_cmplx(cosf(t),sinf(t));

	    m = sf_cmplx(0.,0.);
	    for (ix=nx-1; ix >= 0; ix--) {
#ifdef SF_HAS_COMPLEX_H
		m = m*c + dd[ix];
#else
		m = sf_cadd(sf_cmul(m,c),dd[ix]);
#endif
	    }

	    t = w*p*x0;
	    c = sf_cmplx(cosf(t),sinf(t));
            
            /* transpose */
#ifdef SF_HAS_COMPLEX_H
	    cm[ip][iw] = m*c; 
#else
	    cm[ip][iw] = sf_cmul(m,c);
#endif 
	}
    }

    for (ip=0; ip < np; ip++) { /* loop over slopes */
	/* FFT to time */
	kiss_fftri(invs,(const kiss_fft_cpx *) cm[ip], tt);
		
	for (it=0; it < nt; it++) {
	    pp[it][ip]=tt[it];
	}
    }
    
    /*** dip correction ***/
    for (it=0; it < nt; it++) { /* loop over depth */
	dq = sqrtf(1+dip[it]);
	for (ip=0; ip < np; ip++) {
	    p = (atanf(tanf(p0+ip*dp)*dq)-p0)/dp;
	    /* linear interpolation */
	    iq = floorf(p);
	    if (iq >= 0 && iq < np-1) {
		p -= iq;
		ang[ip][it]=pp[it][iq]*(1.-p)+pp[it][iq+1]*p;
	    } else {
		ang[ip][it]=0.;
	    }
	}
    }
}
