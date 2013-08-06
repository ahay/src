/* Post-stack 2-D velocity continuation. */
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

#include "fint1.h"
#include <math.h>
#include <rsf.h>

#ifdef SF_HAS_FFTW
#include <fftw3.h>
#endif

int main(int argc, char* argv[])
{
    fint1 str, istr;
    bool verb;
    int i1,i2, n1,n2,n3, nw, nx,nv, ix,iv, next;
    float d1,o1,d2,o2, eps, w,x,k, v0,v2,v,dv, dx, t, x0;
    float *trace=NULL, *strace=NULL;
    sf_complex *ctrace=NULL, **cstack=NULL, shift;
    sf_file in=NULL, out=NULL;

#ifdef SF_HAS_FFTW
    fftwf_plan forw, invs;
#else
    kiss_fftr_cfg forw, invs;
#endif

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    if (!sf_getfloat("eps",&eps)) eps=0.01; /* regularization */
    if (!sf_getint("pad",&n2)) n2=n1; /* padding for stretch */
    if (!sf_getint("pad2",&n3)) n3=2*kiss_fft_next_fast_size((n2+1)/2);
    /* padding for FFT */

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */

    nw = n3/2+1;

    strace = sf_floatalloc(n3);
    ctrace = sf_complexalloc(nw);

#ifdef SF_HAS_FFTW
    forw = fftwf_plan_dft_r2c_1d(n3, strace, (fftwf_complex *) ctrace,
				 FFTW_ESTIMATE);
    invs = fftwf_plan_dft_c2r_1d(n3, (fftwf_complex *) ctrace, strace,
				 FFTW_ESTIMATE);
#else
    forw = kiss_fftr_alloc(n3,0,NULL,NULL);
    invs = kiss_fftr_alloc(n3,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("KISS FFT allocation error");
#endif

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;  
    o2 = o1*o1;

    if(!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    d2 = o1+(n1-1)*d1;
    d2 = (d2*d2 - o2)/(n2-1);

    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    /* velocity steps */
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    /* velocity step size */
    if (!sf_getfloat("v0",&v0) && 
	!sf_histfloat(in,"v0",&v0)) sf_error("Need v0=");
    /*( v0 starting velocity )*/

    if(!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if(!sf_histfloat(in,"o2",&x0)) x0=0.;

    sf_putfloat(out,"o2",v0+dv);
    sf_putfloat(out,"d2",dv);
    sf_putint(out,"n2",nv);

    sf_putstring(out,"label2","Velocity");

    dx *= 2.*SF_PI;
    x0 *= 2.*SF_PI;

    trace = sf_floatalloc(n1);
    cstack = sf_complexalloc2(nw,nv);

    if (!sf_getint("extend",&next)) next=4;
    /* trace extension */
    str = fint1_init(next,n1,0);
    istr = fint1_init(next,n2,0);

    for (ix=0; ix < nx; ix++) {
	x = x0+ix*dx; 
	x *= x;

	if (verb) sf_warning("wavenumber %d of %d;", ix+1,nx);
	k = x * 0.25 * 0.25 * 0.5;

	for (iv=0; iv < nv; iv++) {
	    for (i1=0; i1 < nw; i1++) {
		cstack[iv][i1] = sf_cmplx(0.,0.);
	    }
	}

	sf_floatread(trace,n1,in);
	for (i1=0; i1 < n1; i1++) {
	    trace[i1] /= n1;
	}
	fint1_set(str,trace);

	for (i2=0; i2 < n2; i2++) {
	    t = o2+i2*d2;
	    t = sqrtf(t);
	    t = (t-o1)/d1;
	    i1 = t;
	    if (i1 >= 0 && i1 < n1) {
		strace[i2] = fint1_apply(str,i1,t-i1,false);
	    } else {
		strace[i2] = 0.;
	    }
	}
	
	for (i2=n2; i2 < n3; i2++) {
	    strace[i2] = 0.;
	}

#ifdef SF_HAS_FFTW
	fftwf_execute(forw);
#else
	kiss_fftr(forw,strace, (kiss_fft_cpx *) ctrace);
#endif
	ctrace[0]=sf_cmplx(0.,0.); /* dc */

	for (iv=0; iv < nv; iv++) {
	    v = v0 + (iv+1)*dv;
	    v2 = k * ((v0*v0) - (v*v));

	    for (i2=1; i2 < nw; i2++) {
		w = i2*SF_PI/(d2*n3);
		w = v2/w;

		shift = sf_cmplx(cosf(w),sinf(w));

#ifdef SF_HAS_COMPLEX_H
		cstack[iv][i2] += ctrace[i2] * shift;
#else
		cstack[iv][i2] = sf_cadd(cstack[iv][i2],sf_cmul(ctrace[i2],shift));
#endif
	    } /* w */
	} /* v */
 
	for (iv=0; iv < nv; iv++) {
	    for (i2=1; i2 < nw; i2++) {
		w = i2*SF_PI/(d2*n3);
		w *= o2;
		shift = sf_cmplx(cosf(w),sinf(w));

#ifdef SF_HAS_COMPLEX_H
		ctrace[i2] = cstack[iv][i2] * shift;
#else
		ctrace[i2] = sf_cmul(cstack[iv][i2],shift);
#endif
	    }

#ifdef SF_HAS_FFTW
	    fftwf_execute(invs);
#else
	    kiss_fftri(invs,(const kiss_fft_cpx *) ctrace, strace);
#endif
	    fint1_set(istr,strace);
	    
	    for (i1=0; i1 < n1; i1++) {
		t = o1+i1*d1;
		t = t*t;
		t = (t-o2)/d2;
		i2 = t;
		if (i2 >= 0 && i2 < n2) {
		    trace[i1] = fint1_apply(istr,i2,t-i2,false);
		} else {
		    trace[i1] = 0.;
		}
	    }

	    sf_floatwrite (trace,n1,out);
	} /* v 2 */
    } /* x */
    if (verb) sf_warning(".");

    exit (0);
}
