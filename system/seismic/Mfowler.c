/* 2-D velocity-domain imaging (Fowler DMO + Stolt migration).

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
#include "warp2.h"

int main(int argc, char* argv[])
{
    int nk, nw, nv, ik, iw, iv, n;
    float dk, dw, dv, k0, w0, v0, k, w, v, p, vm, wm, s, eps;
    float **stak, **wstr, **vstr, **migr;
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

    stak = sf_floatalloc2(nv,nw);
    wstr = sf_floatalloc2(nv,nw);
    vstr = sf_floatalloc2(nv,nw);
    migr = sf_floatalloc2(nv,nw);

    n = nw*nv;

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    warp2_init(nv, v0, dv,
	       nw, w0, dw,
	       nv, nw, eps);

    for (ik=0; ik < nk; ik++) {
	k = k0+ik*dk;
	sf_floatread(stak[0],n,inp);

	for(iw=0; iw < nw; iw++) {
	    w = w0+iw*dw;
	    if (fabsf(w) < dw) w=SF_SIG(w)*dw;
		
	    p = 0.5*k/w;

	    for (iv=0; iv < nv; iv++) {
		v = v0+iv*dv;
		
		vm = v/hypotf(1.0f,v*p);
		s = vm*p;
		
		if (s < 1.0f) {
		    wm = w*sqrtf(1.0f-s*s);
		    wstr[iw][iv] = wm;
		} else {
		    wstr[iw][iv] = w0-10*dw;
		}

		vstr[iw][iv] = vm;
	    }
	}

	warp2(stak,vstr,wstr,migr);	
	sf_floatwrite(migr[0],n,out);
    }

    exit(0);
}
