/* Offset continuation by finite differences 
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

#include <rsf.h>

#include "ctridiagonal.h"

int main(int argc, char* argv[])
{
    int nw,nh,nx, iw,ix,ih, k;
    float dw, h0,dh,dx, w,w2, h,h2,dh2,hdh;
    float complex diag, diag2, *in, *out, offd, offd2;
    ctris slv;
    sf_file input, output;

    sf_init (argc,argv);
    input = sf_input("in");
    output = sf_output("out");

    if (SF_COMPLEX != sf_gettype(input)) sf_error("Need complex input");

    if (!sf_histint(input,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histint(input,"n2",&nw)) sf_error("No n2= in input");
    if (!sf_histfloat(input,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(input,"d2",&dw)) sf_error("No d2= in input");
    dw *= 2.*SF_PI;

    if (!sf_getint("nh",&nh)) sf_error("Need nh=");
    /* Number of steps in offset */
    if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
    /* Offset step size */
    if (!sf_getfloat("h0",&h0)) sf_error("Need h0=");
    /* Initial offset */

    dh /= dx;
    h0 /= dx;
    dh2 = dh*dh;

    in = sf_complexalloc(nx);
    out = sf_complexalloc(nx);

    slv = ctridiagonal_init (nx);

    for (iw=0; iw < nw; iw++) {
	sf_warning("frequency %d of %d",iw+1,nw);

	w = iw*dw; 
	w2 = w*w;

	sf_complexread(out,nx,input);
	for (ih=0; ih < nh; ih++) {
	    for (ix=0; ix < nx; ix++) {
		in[ix] = out[ix];
	    }

	    h = h0 + ih*dh; 
	    h2 = h*h; 
	    hdh = 2.*dh*h;

	    diag  = 2.*(w*(12.*(h2 - dh2 - hdh) + 5.*w2) + 
			I*(-27.*(dh2 + hdh + 4.*h2) - 
			   3.*(5. + dh2 + hdh)*w2));       
	    diag2 = 2.*(w*(12.*(h2 + 2.*dh2 + 2.*hdh) + 5.*w2) +
			I*(-27.*(3.*dh2 + 3.*hdh + 4.*h2) + 
			   3.*(-5. + dh2 + hdh)*w2));           
	    offd  = w*(12.*(dh2 + hdh - h2) + w2) +
		I*(27.*(dh2 + hdh + 4.*h2) + 3.*(-1. + dh2 + hdh)*w2);
	    offd2 = w*(-12.*(2.*dh2 + 2.*hdh + h2) + w2) + 
		I*(27.*(3.*dh2 + 3.*hdh + 4.*h2) - 3.*(1. + dh2 + hdh)*w2);

	    ctridiagonal_const_define (slv, diag2, offd2);

	    out[0] = diag*in[0] + offd*in[1];
	    for (k = 1; k < nx - 1; k++) {
		out[k] = diag*in[k] + offd*(in[k+1]+in[k-1]);
	    }
	    out[nx-1] = diag*in[nx-1] + offd*in[nx-2];

	    ctridiagonal_solve (slv, out);
	}
	sf_complexwrite (out,nx,output);
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mfincon.c,v 1.5 2004/06/23 23:31:42 fomels Exp $	 */
