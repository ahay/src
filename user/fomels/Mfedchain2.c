/* Fast explicit diffusion as a chain (2-D) */
/*
  Copyright (C) 2022 University of Texas at Austin
  
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

#include "fedchain2.h"

int main(int argc, char* argv[])
{
    bool verb;
    int i, ic, m1, m2, n, nc, n2, iter, niter, liter, k, rect;
    float step, rsum, rsum2;
    float *xn, *x1, *y1, *dx, *r, *x0, *p;
    sf_file inp, out, wht, smo;

    sf_init(argc,argv);
    inp = sf_input("in");
    smo = sf_input("smooth");
    wht = sf_output("w");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&m1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2",&m2)) sf_error("No n2= in input");

    if (!sf_getint("nc",&nc)) nc=1; /* number of components */

    if (!sf_getint("rect",&rect)) rect=1;
    /* smoothing radius */

    sf_smooth1_init(m1,m2,1,rect,rect);

    n = m1*m2;
    n2 = nc*n;
    xn = sf_floatalloc(n2);
    dx = sf_floatalloc(n2);
    x0 = sf_floatalloc(n2);
    r = sf_floatalloc(n2);
    p = sf_floatalloc(n2);
    
    x1 = sf_floatalloc(n);
    sf_floatread(x1,n,inp);

    y1 = sf_floatalloc(n);
    sf_floatread(y1,n,smo);

    fedchain2_init(m1,m2,nc,x1,xn,xn+n);

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=50;
    /* number of linear iterations */

    for (i=0; i < n; i++) {
	xn[i]=1.0f;
    }
    for (ic=1; ic < nc; ic++) {
	for (i=0; i < n; i++) {
	    xn[ic*n+i] = 0.0f;
	}
    }

    sf_conjgrad_init(n2, n2, n2, n2, 1., 1.e-6, verb, false);

    for (iter=0; iter < niter; iter++) {
	fedchain2_apply(y1,r);

	rsum = 0.0f;
	for (i=0; i < n2; i++) {
	    r[i] = -r[i];
	    rsum += r[i]*r[i];
	}

	/*
	sf_solver(fedchain2_lop,sf_cgstep,n2,n2,dx,r,liter,"verb",verb,"end");
	sf_cgstep_close(); */

	sf_conjgrad(NULL, fedchain2_lop, sf_smooth1_lop, p, dx, r, liter);
	
	for (i=0; i < n2; i++) {
	    x0[i] = xn[i];
	}

	/* line search */
	step = 1.0f;
	for (k=0; k < 8; k++) {
	    for (i=0; i < n2; i++) {
		xn[i] = x0[i] + step*dx[i];
	    }
		
	    fedchain2_apply(y1,r);		
	    rsum2 = 0.0f;
	    for (i=0; i < n2; i++) {
		r[i] = -r[i];
		rsum2 += r[i]*r[i];
	    }
		
	    if (rsum2 < rsum) break;
	    
	    step *= 0.5;
	}
    }
    
    sf_floatwrite(xn,n,wht);
    fedchain2(y1);
    sf_floatwrite(y1,n,out);

    exit(0);
}
