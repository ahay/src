/* 2-D velocity-domain DMO (Fowler DMO).

Input: 2-D cosft of constant-velocity stacks (v,w,k).
*/
/*
  Copyright (C) 2013 University of Texas at Austin

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
    sf_map4 mo;
    int nk, nw, nv, ik, iw, iv;
    float dk, dw, dv, k0, w0, v0, k, w, v, p, eps;
    float *stack, *stack2, *vstr;
    sf_file inp, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");

    if (!sf_histint(inp,"n1",&nv)) sf_error("No n1= in input");
    if (!sf_histfloat(inp,"d1",&dv)) sf_error("No d1= in input");
    if (!sf_histfloat(inp,"o1",&v0)) sf_error("No o1= in input");
    
    if (!sf_histint(inp,"n2",&nw)) sf_error("No n2= in input");
    if (!sf_histfloat(inp,"d2",&dw)) sf_error("No d2= in input");
    if (!sf_histfloat(inp,"o2",&w0)) sf_error("No o2= in input");

    if (!sf_histint(inp,"n3",&nk)) sf_error("No n3= in input");
    if (!sf_histfloat(inp,"d3",&dk)) sf_error("No d3= in input");
    if (!sf_histfloat(inp,"o3",&k0)) sf_error("No o3= in input");

    stack  = sf_floatalloc(nv);
    stack2 = sf_floatalloc(nv);
    vstr = sf_floatalloc(nv);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    mo = sf_stretch4_init (nv, v0, dv, nv, eps);

    for (ik=0; ik < nk; ik++) {
	k = k0+ik*dk;

	for(iw=0; iw < nw; iw++) {
	    sf_floatread(stack,nv,inp);

	    w = w0+iw*dw;
	    if (fabsf(w) < dw) w=SF_SIG(w)*dw;
		
	    p = 0.5*k/w;

	    for (iv=0; iv < nv; iv++) {
		v = v0+iv*dv;
		
		vstr[iv] = v/hypotf(1.0f,v*p);
	    }	

	    sf_stretch4_define (mo,vstr,false);
	    sf_stretch4_apply (false,mo,stack,stack2);

	    sf_floatwrite(stack2,nv,out);
	}
    }

    exit(0);
}
