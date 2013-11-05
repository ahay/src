/* Post-stack 2-D velocity continuation in the time-stretch frequency domain. */
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

int main(int argc, char* argv[])
{
    bool verb;
    int nw,nx,nv, iw,ix,iv;
    float w,x,k, v0,v2,v,dv, dx,dw, x0,w0;
    sf_complex *ctrace, shift;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nw)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    
    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
    
    ctrace  = sf_complexalloc(nw);

    if (!sf_histfloat(in,"o1",&w0)) w0=0.;  
    if (!sf_histfloat(in,"d1",&dw)) sf_error("No d1= in input");

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

    sf_shiftdim(in, out, 2);

    dx *= 2.*SF_PI;
    x0 *= 2.*SF_PI;

    dw *= 2.*SF_PI;
    w0 *= 2.*SF_PI;

    for (ix=0; ix < nx; ix++) {
	x = x0+ix*dx; 
	x *= x;

	if (verb) sf_warning("wavenumber %d of %d;", ix+1,nx);
	k = x * 0.25 * 0.25;

	sf_complexread(ctrace,nw,in);

	for (iv=0; iv < nv; iv++) {
	    v = v0 + (iv+1)*dv;
	    v2 = k * ((v0*v0) - (v*v));

	    for (iw=0; iw < nw; iw++) {
		w = w0+iw*dw;

		if (fabsf(w) < dw) w = SF_SIG(w)*dw;

		w = v2/w;
		
		shift = sf_cmplx(cosf(w),sinf(w));
		
#ifdef SF_HAS_COMPLEX_H
		ctrace[iw] *= shift;
#else
		ctrace[iw] = sf_cmul(ctrace[iw],shift);
#endif
	    } /* w */

	    sf_complexwrite (ctrace,nw,out);
	} /* v  */
    } /* x */
    if (verb) sf_warning(".");

    exit (0);
}
