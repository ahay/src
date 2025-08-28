/* Adaptive time-frequency pseudo transform using streaming algorithm. */
/*
  Copyright (C) 2025 Jilin University
  
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
#include "ltfdiv.h"
#include <string.h>

int main(int argc, char* argv[])
{
    int n1, n2, nt, nw, n1w, i1, i2, iw, rect1, *rects;
    float d1, o1, w0, dw, lam, eps0, tol;
    float *inpf, *outf, **bases, *epss,*freq, **rat;
    bool inv, verb, isabs, ismooth, icenter;
    sf_complex *outp;
    sf_file in, out, freqs, feps, frect;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getstring("freq")) freqs=NULL;
    else freqs = sf_input("freq");
    
    if (!sf_getstring("epss")) feps=NULL;
    else feps = sf_input("epss");
    
    if (!sf_getstring("rects")) frect=NULL;
    else frect = sf_input("rects");

    
    if(NULL!=freqs){
        if(sf_gettype(freqs)!=SF_FLOAT) sf_error("Need float frequencies");
        sf_histint(freqs,"n1",&nw);
        dw=-1;w0=-1;
        freq = sf_floatalloc(nw);
        sf_floatread(freq,nw,freqs);
    }

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1 = 0. ;

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    /* if y, do inverse transform */
    if (!sf_getbool("smooth",&ismooth)) ismooth = false;
    /* verbosity flag */
    if (!sf_getbool("abs",&isabs)) isabs = false;

    if (!sf_getbool("center",&icenter)) icenter = true;
    /* verbosity flag */
    if (!sf_getfloat("lambda",&lam)) lam=1.0f;
    /* regularization trade-off parameter */
    lam *=lam;

    if (!sf_getfloat("tol",&tol)) tol=10.f;
    if (!sf_getfloat("eps",&eps0)) eps0=0.95f;
    /* smoothing parameter */
    if (eps0>1.) eps0 =1.f; 
    eps0 *= eps0;

    if (!sf_getint("rect1",&rect1)) rect1=1;
    /* smooth radius */

    if (!inv) {
        if(NULL ==freqs){
	    if (!sf_getint("nw",&nw)) { /* number of frequencies */
		nt = 2*kiss_fft_next_fast_size((n1+1)/2);
		nw = nt/2+1;
		dw = 1./(nt*d1);
		w0 = 0.;
	    } else {
		if (!sf_getfloat("dw",&dw)) {
		    /* frequency step */
		    nt = 2*kiss_fft_next_fast_size((n1+1)/2);
		    dw = 1./(nt*d1);
		}
		if (!sf_getfloat("w0",&w0)) w0=0.;
		/* first frequency */
	    }
        }
        n2 = sf_leftsize(in,1);
        sf_shiftdim(in, out, 2);
        sf_putint(out,"n2",nw);
        sf_putfloat(out,"d2",dw);
        sf_putfloat(out,"o2",w0);
        sf_putstring(out,"label2","Frequency");
        sf_putstring(out,"unit2","Hz");
        if (isabs)  sf_settype(out,SF_FLOAT);
        else        sf_settype(out,SF_COMPLEX);
	
    } else {
        n2 = sf_leftsize(in,2);
        if(NULL ==freqs){
	    if (!sf_histint(in,"n2",&nw)) sf_error("No n2= in input");
	    if (!sf_histfloat(in,"d2",&dw)) sf_error("No d2= in input");
	    if (!sf_histfloat(in,"o2",&w0)) sf_error("No o2= in input");
        }
        sf_unshiftdim(in, out, 2);
        sf_settype(out,SF_FLOAT);
    }
    
    n1w = n1*nw;
    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;
    
    inpf = sf_floatalloc(n1);
    bases  = sf_floatalloc2(n1,2*nw);
    outp = sf_complexalloc(n1w);
    if(isabs) outf = sf_floatalloc(n1w);
    else outf = NULL;
    epss = sf_floatalloc(nw);
    rects = sf_intalloc(nw);

    if(!feps){
        for (iw=0; iw<nw; iw++) epss[iw] = eps0;
    }
    else sf_floatread(epss,nw,feps);
    
    if(!frect){
        for (iw=0; iw<nw; iw++) rects[iw] = rect1;
    }
    else sf_intread(rects,nw,frect);
    
    if (!inv) {
        rat = sf_floatalloc2(n1,2*nw);
        memset(*rat, 0, sizeof(float)*n1*2*nw);
    } else {
	rat = NULL;
    }
    
    
    /* Sine and cosine bases */
    for (i1=0; i1<n1; i1++) {
	for (iw=0; iw<nw; iw++){
	    if(NULL==freqs){
		bases[iw][i1]    = cos((w0+iw*dw)*(i1*d1+o1));
		bases[iw+nw][i1]  = sin((w0+iw*dw)*(i1*d1+o1));
	    } else{
		bases[iw][i1]    = cos(2.*SF_PI*(freq[iw])*(i1*d1+o1));
		bases[iw+nw][i1]  = sin(2.*SF_PI*(freq[iw])*(i1*d1+o1));
	    }
	}
    }
    if(!inv) ltfdiv_init(n1, 2*nw, rects, epss, lam, tol, icenter, ismooth);
    
    for (i2=0; i2<n2; i2++){ /* Loop over traces */
	if(verb) sf_warning("trace %d of %d;",i2+1,n2);
        if(!inv){
            sf_floatread(inpf, n1, in);
            ltfdiv(inpf, bases, rat);
	    
            if (!isabs){
		for (iw = 0; iw < nw; iw++) {
		    for (i1 = 0; i1 < n1; i1++) {
			outp[iw*n1+i1] = sf_cmplx(rat[iw][i1],
						  rat[iw+nw][i1]);
		    }
		    
		}
		sf_complexwrite(outp,n1w,out);
            }/* not abs */
            else{
		for (iw = 0; iw < nw; iw++) {
		    for (i1 = 0; i1 < n1; i1++) {
			outf[iw*n1+i1] = sqrtf(rat[iw][i1]*rat[iw][i1]+	\
					       rat[iw+nw][i1]*rat[iw+nw][i1]);
		    }
		}
		sf_floatwrite(outf,n1w,out);
            }/* abs */
            memset(*rat, 0, sizeof(float)*n1*2*nw);
        } else {
            // memset(outp, 0, sizeof(sf_complex)*n1w);
            memset(inpf, 0, sizeof(float)*n1);
            sf_complexread(outp,n1w,in);
            for (iw=0; iw < nw; iw++) {
		for (i1=0; i1 < n1; i1++) {
		    inpf[i1] += crealf(outp[iw*n1+i1])*bases[iw][i1]+
			cimagf(outp[iw*n1+i1])*bases[iw+nw][i1];
		}
            }
            sf_floatwrite(inpf,n1,out);
        }
    }
    sf_warning(".");
    if(!inv) ltfdiv_close();
    exit(0);
}
