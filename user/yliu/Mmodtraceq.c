/* Generate single trace with Q attenuation for viscoelastic media. */
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
    bool  verb;
    int   n1, nc, fm, *q=NULL, *at=NULL, nw, nt, wr, i1, iw, i, nwc;
    float d1, o1, dw, w0, w;
    float *R, *trace, *T, *tra;
    /* char  key[7]; */
    kiss_fft_cpx *bs, *ibs, *enat, *eenat;
    kiss_fftr_cfg icfg;
    
    sf_file out;

    sf_init(argc, argv);

    out = sf_output("out");

    /* verbosity */
    if (!sf_getbool("verb",&verb)) verb = false;

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
	q = sf_intalloc(nc);

	/* snprintf(key,6,"at"); */
	if (!sf_getints("at",at,nc)) {
	    /* (at=[at1,at2,...] layer quality factor) */
	    sf_error("need at");
	}

	/* snprintf(key,6,"q"); */
	if (!sf_getints("q",q,nc)) {
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
    wr = 500.;
    /* qr = 500.; */

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"o1",o1);

    nwc=nw*nc;

    R = sf_floatalloc(nw);
    enat = (kiss_fft_cpx*)sf_complexalloc(nwc);
    eenat = (kiss_fft_cpx*)sf_complexalloc(nwc);
    bs = (kiss_fft_cpx*)sf_complexalloc(nwc);
    ibs = (kiss_fft_cpx*)sf_complexalloc(nw);
    tra = sf_floatalloc(nt);
    trace = sf_floatalloc(n1);
    T = sf_floatalloc(nc);
    
    if (!sf_getint("fm",&fm)) fm = 50;
    /* dominant frequency of Ricker wavelet */

    /* Ricker wavelet(initial wavelet) */
    for (iw = 0; iw < nw; iw++) {
	w = w0 + iw * dw;
	R[iw] = (2/sqrt(SF_PI))*((w*w)/(fm*fm))*exp(-((w*w)/(fm*fm)));
    }

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
	    bs[i*nw+iw] = sf_crmul(enat[i*nw+iw],R[iw]);
	}
    }

    icfg = kiss_fftr_alloc(nt,1,NULL,NULL);

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

    kiss_fftri(icfg,ibs,tra);

    for (i1 = 0; i1 < n1; i1++) {
      trace[i1] = tra[i1];
    }

    sf_floatwrite(trace,n1,out);
    
}
