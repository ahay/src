/* Fitting a polynomial by least-squares. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
    int n1, n2, i2, i1, ic, id, nc, n;
    float *dat, *crd, **func, **mat, *rhs, *sol, *dreg, x, xp, eps, o, d;
    sf_file inp, coord, coef, out, reg;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_getfloat("eps",&eps)) eps=0.0f;
    /* regularization parameter */
 
    if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(inp,1);

    if (!sf_getint("np",&nc)) nc=1; 
    /* polynomial order */
    nc++;

    if (NULL != sf_getstring("coef")) {
        /* (optional) coefficients */
	coef = sf_output("coef"); 
	sf_putint(coef,"n1",nc);
    } else {
	coef = NULL;
    }


    dat = sf_floatalloc(n1);
    crd = sf_floatalloc(n1);
    func = sf_floatalloc2(n1,nc);

    if (NULL != sf_getstring("coord")) {
	/* coordinates */
	coord = sf_input("coord");
    } else {
	coord = NULL;
	if (!sf_histfloat(inp,"d1",&d)) d=1.0f;
	if (!sf_histfloat(inp,"o1",&o)) o=0.0f;
	for (i1=0; i1 < n1; i1++) {
	    crd[i1] = o+i1*d;
	}
    }

    if (NULL != sf_getstring("reg")) {
	/* (optional) regularly sampled */
	reg = sf_output("reg");
	if (!sf_getint("n1",&n)) sf_error("Need n1=");
	/* number of samples for regularly sampled */
	if (!sf_getfloat("d1",&d)) sf_error("Need d1=");
	/* sampling for regularly sampled */
	if (!sf_getfloat("o1",&o)) sf_error("Need o1=");
	/* origin for regularly sampled */
	sf_putint(reg,"n1",n);
	sf_putfloat(reg,"d1",d);
	sf_putfloat(reg,"o1",o);

	dreg = sf_floatalloc(n);
    } else {
	reg = NULL;
	dreg = NULL;
    }

    sf_gaussel_init(nc);
    sol = sf_floatalloc(nc);
    rhs = sf_floatalloc(nc);
    mat = sf_floatalloc2(nc,nc);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(dat,n1,inp);
	if (NULL != coord) sf_floatread(crd,n1,coord);

	for (i1=0; i1 < n1; i1++) {
	    x = crd[i1];
	    xp = 1.0f;
	    for (ic=0; ic < nc; ic++) {
		func[ic][i1] = xp;
		xp *= x;
	    }
	}
	
	/* compute A'A matrix */
	for (ic=0; ic < nc; ic++) {
	    for (id=0; id <= ic; id++) {
		mat[ic][id] = 0.;
		for (i1=0; i1 < n1; i1++) {
		    mat[ic][id] += func[ic][i1]*func[id][i1];
		}
		mat[id][ic] = mat[ic][id];
	    }
	}

	if (eps > 0.0f) { 
            /* regularization: use A'A + eps*I */
	    for (ic=0; ic < nc; ic++) {
		mat[ic][ic] += eps;
	    }
	}

	/* compute A'd */
	for (ic=0; ic < nc; ic++) {
	    rhs[ic] = 0.;
	    for (i1=0; i1 < n1; i1++) {
		rhs[ic] += func[ic][i1]*dat[i1];
	    }
	}

	/* inversion */
	sf_gaussel_solve(mat,rhs,sol);

	if (NULL != coef) sf_floatwrite(sol,nc,coef);

	/* compute Ac */
	for (i1=0; i1 < n1; i1++) {
	    dat[i1] = 0.;
	    for (ic=0; ic < nc; ic++) {
		dat[i1] += func[ic][i1]*sol[ic];
	    }
	}

	sf_floatwrite(dat,n1,out);

	if (NULL != reg) {
	    for (i1=0; i1 < n; i1++) {
		dreg[i1] = 0.;
		x = o+i1*d;
		xp = 1.0f;
		for (ic=0; ic < nc; ic++) {
		    dat[i1] += xp*sol[ic];
		    xp *= x;
		}
	    }
	    sf_floatwrite(dreg,n,reg);
	}
    }

    exit(0);
}
