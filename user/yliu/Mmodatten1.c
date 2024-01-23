/* 1D attenuation modeling according to modified Kolsky model. */
/*
  Copyright (C) 2022 Jilin University
  
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
#include <math.h>

int main(int argc, char *argv[])
{
    int   n1, nc, fm, *at=NULL;
    float *q=NULL, mag;
    int   nw, nt, wr, i1, iw, i, nwc;
    float d1, o1, dw, w0, w;
    float *trace, *T;
    float *min, *minp;
    kiss_fft_cpx *R, *bs, *ibs, *enat, *eenat;
    kiss_fftr_cfg cfg, icfg;   
    sf_file out;

    sf_init(argc, argv);

    out = sf_output("out");

    /* basic parameters */
    if (!sf_getint("n1",&n1)) n1 = 1000;
    /* size of time axis */  
    
    if (!sf_getfloat("d1",&d1)) d1 = 0.001;
    /* sampling on time axis */
    
    if (!sf_getfloat("o1",&o1)) o1 = 0.;
    /* origin on time axis */
    
    if (!sf_getint("nc",&nc)) nc = 1;
    /* number of layer */
    if (nc < 1) sf_error("Need nq >= 1");

    if (nc >= 1) {
	at = sf_intalloc(nc);
	q = sf_floatalloc(nc);

	/* snprintf(key,6,"at"); */
	if (!sf_getints("at",at,nc)) {
	    /* (at=[at1,at2,...] layer quality factor) */
	    sf_error("need at");
	}

	/* snprintf(key,6,"q"); */
	if (!sf_getfloats("q",q,nc)) {
	    /* (q=[q1,q2,...] layer quality factor) */
	    for (i = 0; i < nc; i++) {
		q[i] = 100;
	    }
	}
    }
    
    nt = 2*kiss_fft_next_fast_size((n1+1)/2);
    nw = nt/2+1;
    dw = 1./(nt*d1);
    w0 = 0.;
    wr = 500;

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"o1",o1);

    nwc=nw*nc;

    R = (kiss_fft_cpx*)sf_complexalloc(nw);
    min = sf_floatalloc(n1);
    minp = sf_floatalloc(n1);
    enat = (kiss_fft_cpx*)sf_complexalloc(nwc);
    eenat = (kiss_fft_cpx*)sf_complexalloc(nwc);
    bs = (kiss_fft_cpx*)sf_complexalloc(nwc);
    ibs = (kiss_fft_cpx*)sf_complexalloc(nw);
    trace = sf_floatalloc(n1);
    T = sf_floatalloc(nc);
    
    if (!sf_getint("fm",&fm)) fm = 40;
    /* dominant frequency of Ricker wavelet */

    if (!sf_getfloat("mag",&mag)) mag = 1.;

    /* min-phase */
    for (i1=0;i1<n1;i1++) {
    	min[i1]=2*SF_PI*fm*i1*d1;
    	minp[i1]=mag*exp(-min[i1]*min[i1]/25)*sinf(min[i1]);
    }
    cfg = kiss_fftr_alloc(nt,0,NULL,NULL);
    kiss_fftr(cfg,minp,R);

    /* seismis wavelet attenuation */

    /* frequency dispersion */
    for (i = 0; i < nc; i++ ) {
    	T[i] = at[i] * d1;
    }
    for (i = 0; i < nc; i++ ) {
    	for (iw = 0; iw < nw; iw++) {
    	    w = w0 + iw * dw;
    	    eenat[i*nw+iw].r = 0;
    	    eenat[i*nw+iw].i = 0;
    	    if (w == 0) {
    		enat[i*nw+iw].r = 1;
    		enat[i*nw+iw].i = 0;
    	    } else {
    		if (i == 0) {
    		    enat[i*nw+iw].r = exp(-SF_PI*w*T[i]/q[i])
    		    	*cos(2*SF_PI*w*T[i]/(pow(w/wr,1/(SF_PI*q[i]))));
    		    enat[i*nw+iw].i = exp(-SF_PI*w*T[i]/q[i])
    		    	*sin(-2*SF_PI*w*T[i]/(pow(w/wr,1/(SF_PI*q[i]))));
    		} else {
    		    eenat[i*nw+iw].r = exp(-SF_PI*w*(T[i]-T[i-1])/q[i])
    		    	*cos(2*SF_PI*w*(T[i]-T[i-1])/(pow(w/wr,1/(SF_PI*q[i]))));
    		    eenat[i*nw+iw].i = exp(-SF_PI*w*(T[i]-T[i-1])/q[i])
    		    	*sin(-2*SF_PI*w*(T[i]-T[i-1])/(pow(w/wr,1/(SF_PI*q[i]))));
    		    enat[i*nw+iw] = sf_cmul(enat[(i-1)*nw+iw],eenat[i*nw+iw]);
    		}
    	    }
    	}
    }

    for (i = 0; i < nc; i++ ) {
	for (iw = 0; iw < nw; iw++) {
	    bs[i*nw+iw] = sf_cmul(enat[i*nw+iw],R[iw]);
	}
    }

    icfg = kiss_fftr_alloc(nt,1,NULL,NULL);

    for (i1 = 0; i1 < n1; i1++) {
	trace[i1] = 0.;
    }
    
    for (iw = 0; iw < nw; iw++) {
    	ibs[iw].r=0.;
    	ibs[iw].i=0.;
    }
    
    for (i = 0; i < nc; i++ ) {
    	for (iw = 0; iw < nw; iw++) {
    	    ibs[iw].r = ibs[iw].r + bs[i*nw+iw].r;
    	    ibs[iw].i = ibs[iw].i + bs[i*nw+iw].i;
    	}
    }

    kiss_fftri(icfg,ibs,trace);
    
    for (i1 = 0; i1 < n1; i1++) {
      trace[i1] = trace[i1]/nt;
    }
    
    sf_floatwrite(trace,n1,out);
    
}
