/* Velocity continuation.
*/
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

int main(int argc, char* argv[])
{
    fint1 str, istr;
    int i1,i2, n1,n2,n3, nw, nx,ny,nv,nh, ix,iy,iv,ih, next;
    float d1,o1,d2,o2, eps, w,x,y,k, v0,v2,v,v1,dv, dx,dy, h0,dh,h, t;
    float *trace, *strace;
    float complex *ctrace, **cstack;
    kiss_fftr_cfg forw, invs;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nx)) sf_error("No n3= in input");
    if (!sf_histint(in,"n4",&ny)) ny=1;

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    if (!sf_getint("pad",&n2)) n2=n1;
    if (!sf_getint("pad2",&n3)) n3=n2;

    if (n3%2) n3++;
    nw = n3/2+1;
    forw = kiss_fftr_alloc(n3,0,NULL,NULL);
    invs = kiss_fftr_alloc(n3,1,NULL,NULL);
    if (NULL == forw || NULL == invs) 
	sf_error("KISS FFT allocation error");

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;  
    o2 = o1*o1;
 
    if(!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    d2 = o1+(n1-1)*d1;
    d2 = (d2*d2 - o2)/(n2-1);
    
    if (!sf_getint("nv",&nv)) sf_error("Need nv=");
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");
    if (!sf_getfloat("v0",&v0) && 
	!sf_histfloat(in,"v0",&v0)) sf_error("Need v0=");

    if(!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input");
    if(!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if(!sf_histfloat(in,"d3",&dx)) sf_error("No d3= in input");
    if(!sf_histfloat(in,"d4",&dy)) dy=dx;

    sf_putfloat(out,"o2",v0+dv);
    sf_putfloat(out,"d2",dv);
    sf_putint(out,"n2",nv);

    sf_putstring(out,"label2","Velocity (km/s)");

    dx *= 2.*SF_PI; 
    dy *= 2.*SF_PI; 

    trace = sf_floatalloc(n1);
    strace = sf_floatalloc(n3);
    ctrace = sf_complexalloc(nw);
    cstack = sf_complexalloc2(nw,nv);

    if (!sf_getint("extend",&next)) next=4;
    /* trace extension */
    str = fint1_init(next,n1,0);
    istr = fint1_init(next,n2,0);

    for (iy=0; iy < ny; iy++) {
	y = iy*dy; 
	y *= y;
	for (ix=0; ix < nx; ix++) {
	    x = ix*dx; 
	    x *= x;

	    sf_warning("wavenumber %d of %d and %d of %d", ix+1,nx, iy+1,ny);
	    k = (x+y) * 0.25 * 0.25 * 0.5;

	    for (iv=0; iv < nv; iv++) {
		for (i1=0; i1 < nw; i1++) {
		    cstack[iv][i1] = 0.;
		}
	    }

	    for (ih=0; ih < nh; ih++) {
		h = h0 + ih*dh;
		h = 4. * h * h;

		sf_floatread(trace,n1,in);
		for (i1=0; i1 < n1; i1++) {
		    trace[i1] /= (n1*nh);
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
		
		kiss_fftr(forw,strace, (kiss_fft_cpx *) ctrace);       	
		ctrace[0]=0.; /* dc */

		for (iv=0; iv < nv; iv++) {
		    v = v0 + (iv+1)*dv;
		    v1 = h * (1./(v*v) - 1./(v0*v0)) * 2.;
		    v2 = k * ((v0*v0) - (v*v));

		    for (i2=1; i2 < nw; i2++) {
			w = i2*SF_PI/(d2*n3);
 
			cstack[iv][i2] += 
			    ctrace[i2] * cexpf(I*(v2/w+(v1-o2)*w));
		    } /* w */
		} /* v */
	    } /* h */
 
	    for (iv=0; iv < nv; iv++) {
		for (i2=1; i2 < nw; i2++) {
		    w = i2*SF_PI/(d2*n3);
		    ctrace[i2] = cstack[iv][i2] * cexpf(I*o2*w);
		}

		kiss_fftri(invs,(const kiss_fft_cpx *) ctrace, strace);
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
    } /* y */

    exit (0);
}

/* 	$Id$	 */
