/* Diffraction focusing test. */
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

int main(int argc, char* argv[]) 
{
    int it, nt, ix, nx, i2, n2, n3, iw, nw, next;
    float v0, v1, dt, t0, o2, d2, t, dx, k, w, dw;
    float **data, *trace, *strace;
    sf_complex *ctrace, shift;
    fint1 str, istr;
    kiss_fftr_cfg forw, invs;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&t0)) t0=0.;  

    if (!sf_histint(inp,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d2",&dx)) sf_error("No d2= in input");

    if (!sf_getfloat("v0",&v0)) v0=SF_EPS;
    /* initial velocity */

    if (!sf_getfloat("v",&v1)) sf_error("Need v=");
    /* final velocity */

    if (!sf_getint("pad",&n2)) n2=nt; /* padding for stretch */
    if (!sf_getint("pad2",&n3)) n3=2*kiss_fft_next_fast_size((n2+1)/2);
    /* padding for FFT */

    o2 = t0*t0;
    d2 = t0+(nt-1)*dt;
    d2 = (d2*d2 - o2)/(n2-1);

    if (!sf_getint("extend",&next)) next=4;
    /* trace extension */
    str = fint1_init(next,nt,0);
    istr = fint1_init(next,n2,0);

    nw = n3/2+1;
    dw = 2*SF_PI/(n3*d2);

    forw = kiss_fftr_alloc(n3,0,NULL,NULL);
    invs = kiss_fftr_alloc(n3,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("KISS FFT allocation error");

    data = sf_floatalloc2(nt,nx);
    strace = sf_floatalloc(n3);
    ctrace = sf_complexalloc(nw);

    sf_floatread(data[0],nt*nx,inp);

    /* cosine FT in space */
    sf_cosft_init(nx);
    dx = SF_PI/(kiss_fft_next_fast_size(nx-1)*dx);

    for (it=0; it < nt; it++) {
	sf_cosft_frw (data[0], it, nt);
    }

    for (ix=0; ix < nx; ix++) {
	/* loop over wavenumbers */
	k = ix*dx;
	k *= k*(v0*v0 - v1*v1);

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

	/* FFT */

	kiss_fftr(forw,strace, (kiss_fft_cpx *) ctrace);

	/* velocity continuation itself */

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

    sf_floatwrite(data[0],nt*nx,out);

    exit(0);
}
