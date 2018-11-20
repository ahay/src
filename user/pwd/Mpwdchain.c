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

#include "pwdchain.h"
#include "smooth1.h"

int main(int argc, char* argv[])
{
    bool verb;
    int i, ic, m1, m2, n, nc, n2, iter, niter, liter, rect1, rect2, it, nt, nw;
    float *xn, *x1, *y1, *dx, *r, *p;
    sf_file inp, out, dip;

    sf_init(argc,argv);
    inp = sf_input("in");
    dip = sf_output("dip");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    if (!sf_histint(inp,"n1",&m1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n1",&m2)) sf_error("No n2= in input");
    n = m1*m2;
    nt = sf_leftsize(inp,2);

    if (!sf_getint("nc",&nc)) nc=1; /* number of components */

    sf_putint(dip,"n3",nc);
    sf_shiftdim(inp, dip, 3);
    
    n2 = (2*nc-1)*n;

    xn = sf_floatalloc(n2);
    dx = sf_floatalloc(n2);
    p  = sf_floatalloc(n2);
    
    x1 = sf_floatalloc(n);
    y1 = sf_floatalloc(n);

    if (!sf_getint("order",&nw)) nw=1; /* PWD order */

    pwdchain_init(m1,m2,nw,nc,x1,xn,xn+n*nc);
    r = sf_floatalloc(n*nc);

    if (!sf_getbool("verb",&verb)) verb=(bool) (1 == nt);
    /* verbosity flag */
    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=50;
    /* number of linear iterations */

    if (!sf_getint("rect1",&rect1)) rect1=1;
    if (!sf_getint("rect2",&rect2)) rect2=1;
    /* smoothing radius */

    smooth1_init(m1,m2,nc,rect1,rect2);

    for (it=0; it < nt; it++) {
	sf_warning("slice %d of %d;",it+1,nt);
	sf_floatread(x1,n,inp);
	
	/* initialize */
	for (i=0; i < n; i++) {
	    y1[i] = 0.0f;
	}
	

	if (1==nc) {
	    for (i=0; i < n; i++) {
		xn[i] = 0.0f;
	    }
	} else {
	    for (ic=0; ic < nc; ic++) {
		for (i=0; i < n; i++) {
		    xn[ic*n+i] =  -nw + 2*nw*ic/(nc-1);
		    /* distribute from -nw to nw */
		}
	    }
	}
	for (ic=0; ic < nc-1; ic++) {
	    for (i=0; i < n; i++) {
		xn[(nc+ic)*n+i] = 0.0f; 
	    }
	}
	
	sf_cconjgrad_init(n2, n2, n*nc, n*nc, 1., 1.e-6, verb, false);

	for (iter=0; iter < niter; iter++) {
	    pwdchain_apply(y1,r);
	    
	    for (i=0; i < n*nc; i++) {
		r[i] = -r[i];
	    }

	    sf_conjgrad(NULL, pwdchain_lop, smooth1_lop,p,dx,r,liter);

	    for (i=0; i < n2; i++) {
		xn[i] += dx[i];
	    }
	}

	sf_cconjgrad_close();

	sf_floatwrite(xn,n*nc,dip);

	pwdchain(y1);

	sf_floatwrite(y1,n,out);
    }
    sf_warning(".");

    exit(0);
}
