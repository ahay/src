/* Nonstationary Prony by chain of PEFs */
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

#include "pchain.h"
#include "csmooth1.h"

int main(int argc, char* argv[])
{
    bool verb;
    int i, ic, n, nc, n2, iter, niter, liter, rect, it, nt;
    float f;
    sf_complex *xn, *x1, *y1, *dx, *r, *p;
    sf_file inp, out, pef;

    sf_init(argc,argv);
    inp = sf_input("in");
    pef = sf_output("pef");
    out = sf_output("out");

    if (SF_COMPLEX != sf_gettype(inp)) sf_error("Need complex input");
    if (!sf_histint(inp,"n1",&n)) sf_error("No n1= in input");
    nt = sf_leftsize(inp,1);

    if (!sf_getint("nc",&nc)) nc=1; /* number of components */

    sf_putint(pef,"n2",nc);
    sf_shiftdim(inp, pef, 2);
    
    n2 = (2*nc-1)*n;

    xn = sf_complexalloc(n2);
    dx = sf_complexalloc(n2);
    p  = sf_complexalloc(n2);
    
    x1 = sf_complexalloc(n);

    y1 = sf_complexalloc(n);

    pchain_init(n,nc,x1,xn,xn+n*nc);
    r = sf_complexalloc(n*nc);

    if (!sf_getbool("verb",&verb)) verb=(bool) (1 == nt);
    /* verbosity flag */
    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=50;
    /* number of linear iterations */

    if (!sf_getint("rect",&rect)) rect=1;
    /* smoothing in time */

    csmooth1_init(n,nc,rect);

    for (it=0; it < nt; it++) {
	sf_warning("trace %d of %d;",it+1,nt);
	sf_complexread(x1,n,inp);
	
	/* initialize */
	for (i=0; i < n; i++) {
	    y1[i] = sf_cmplx(0.0f,0.0f);
	}
	
	for (ic=0; ic < nc; ic++) {
	    f = SF_PI*ic/nc;
	    for (i=0; i < n; i++) {
		xn[ic*n+i] =  cexpf(sf_cmplx(0.0f,f)); 
		/* distribute around the unit circle */
	    }
	}
	for (ic=0; ic < nc-1; ic++) {
	    for (i=0; i < n; i++) {
		xn[(nc+ic)*n+i] = sf_cmplx(0.0f,0.0f); 
	    }
	}

	sf_cconjgrad_init(n2, n2, n*nc, n*nc, 1., 1.e-6, verb, false);

	for (iter=0; iter < niter; iter++) {
	    pchain_apply(y1,r);
	    
	    for (i=0; i < n*nc; i++) {
		r[i] = -r[i];
	    }

	    sf_cconjgrad(NULL, pchain_lop, csmooth1_lop,p,dx,r,liter);

	    for (i=0; i < n2; i++) {
		xn[i] += dx[i];
	    }
	}

	sf_cconjgrad_close();

	sf_complexwrite(xn,n*nc,pef);

	pchain(y1);

	sf_complexwrite(y1,n,out);
    }
    sf_warning(".");

    exit(0);
}
