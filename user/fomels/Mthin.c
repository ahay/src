/* Sparse deconvolution. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    int i1, n1, nt, i2, n2, n12, nclip, iw, nw, iter, niter;
    float *wvl, *iwvl, eps, w, pclip, clip, **s, **r, *rr, **ar, *err;
    sf_complex **fs, *fr;
    kiss_fftr_cfg cfg, icfg;
    sf_file seis, refl, wave, conv;

    sf_init(argc,argv);
    seis = sf_input("in");
    refl = sf_output("out");
    wave = sf_input("wave");

    if (!sf_histint(seis,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(seis,1);
    n12 = n1*n2;

    if (!sf_histint(wave,"n1",&nw)) sf_error("No n1= in wave");
    nt = 2*(nw-1);

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */

    if (!sf_getfloat("eps",&eps)) eps=0.0001;
    /* regularization for Wiener deconvolution */

    if (!sf_getfloat("pclip",&pclip)) pclip=4;
    /* percentage to threshold */
    nclip = 0.5+n12*(1.-0.01*pclip);
    if (nclip < 0) nclip=0;
    if (nclip >= n12) nclip=n12-1;

    if (NULL != sf_getstring("conv")) {
	conv = sf_output("conv");
	sf_putint(conv,"n1",niter);
	sf_putfloat(conv,"o1",1.0);
	sf_putfloat(conv,"d1",1.0);
	sf_putstring(conv,"label1","Iteration");
	sf_putstring(conv,"unit1","");
	sf_putint(conv,"n2",1);
	sf_putint(conv,"n3",1);
	err = sf_floatalloc(niter);
    } else {
	conv = NULL;
	err = NULL;
    }
    
    wvl = sf_floatalloc(nw);
    iwvl = sf_floatalloc(nw);

    sf_floatread(wvl,nw,wave);

    /* Wiener deconvolution */
    for (iw=0; iw < nw; iw++) {
	w = wvl[iw];
	iwvl[iw] = w/(w*w+eps);
    }

    r = sf_floatalloc2(nt,n2);
    ar = sf_floatalloc2(n1,n2);
    rr = sf_floatalloc(nt);
    s = sf_floatalloc2(nt,n2);

    /* initialize reflectivity with zero */
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < nt; i1++) {
	    r[i2][i1] = 0.;
	}
	sf_floatread(s[i2],n1,seis);
	for (i1=n1; i1 < nt; i1++) {
	    s[i2][i1] = 0.;
	}
    }

    /* Fourier transform */
    fs = sf_complexalloc2(nw,n2);
    fr = sf_complexalloc(nw);

    cfg  = kiss_fftr_alloc(nt,0,NULL,NULL);
    icfg = kiss_fftr_alloc(nt,1,NULL,NULL);

    for (i2=0; i2 < n2; i2++) {
	kiss_fftr (cfg,s[i2],(kiss_fft_cpx *) fs[i2]);
    }

    for (iter=0; iter < niter; iter++) {
	/* res = F d - W F r */ 
	/* r_next = T [ F^{-1} W/(W^2+eps) res + r ] */

	if (NULL != err) err[iter] = 0.;

	for (i2=0; i2 < n2; i2++) {
	    kiss_fftr (cfg,r[i2],(kiss_fft_cpx *) fr);
	    for (iw=0; iw < nw; iw++) {
#ifdef SF_HAS_COMPLEX_H
		fr[iw] = fs[i2][iw]-wvl[iw]*fr[iw];
		if (NULL != err) err[iter] += cabsf(fr[iw]*conjf(fr[iw]));
		fr[iw] *= iwvl[iw];
#else
		fr[iw] = sf_cadd(fs[i2][iw],sf_crmul(fr[iw],-wvl[iw]));
		if (NULL != err) err[iter] += cabsf(sf_cmul(fr[iw],conjf(fr[iw])));
		fr[iw] = sf_crmul(fr[iw],iwvl[iw]);
#endif	       
	    }
	    kiss_fftri(icfg,(kiss_fft_cpx *) fr,rr);
	    for (i1=0; i1 < n1; i1++) {
		r[i2][i1] += rr[i1]/nt;
		ar[i2][i1] = fabsf(r[i2][i1]);
	    }
	}
	 
	if (NULL != err && iter > 0) 
	    sf_warning("iter=%d err=%g",iter,err[iter]/err[0]);

	/* later change thresholding to sharpening */
	clip = sf_quantile(nclip,n12,ar[0]);
	
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		w = r[i2][i1];
		if (w < -clip) {
		    r[i2][i1] = w+clip;
		} else if (w > clip) {
		    r[i2][i1] = w-clip;
		} else {
		    r[i2][i1] = 0.;
		}
	    }
	}
    } /* iter */
    
    if (NULL != err) sf_floatwrite(err,niter,conv);

    for (i2=0; i2 < n2; i2++) {
	sf_floatwrite(r[i2],n1,refl);
    }

    exit(0);
}
