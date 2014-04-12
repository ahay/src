/* From roots to polynomials */
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
    int n2, i2, i1, ic, id, nc, n;
    double xp;
    float *roots, *sol, *dat, x, o, d;
    sf_file inp, coef, out;

    sf_init(argc,argv);
    inp = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(inp,"n1",&nc))   sf_error("No n1= in input");
    nc++;
    n2 = sf_leftsize(inp,1);

    if (NULL != sf_getstring("coef")) {
        /* (optional) coefficients */
	coef = sf_output("coef"); 
	sf_putint(coef,"n1",nc);
    } else {
	coef = NULL;
    }

    roots = sf_floatalloc(nc-1);

    out = sf_output("out");
    if (!sf_getint("n1",&n)) sf_error("Need n1=");
    /* number of samples */
    if (!sf_getfloat("d1",&d)) sf_error("Need d1=");
    /* sampling */
    if (!sf_getfloat("o1",&o)) sf_error("Need o1=");
    /* origin */
    sf_putint(out,"n1",n);
    sf_putfloat(out,"d1",d);
    sf_putfloat(out,"o1",o);

    dat = sf_floatalloc(n);
    sol = sf_floatalloc(nc);
    sol[0]=1.0;

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(roots,nc-1,inp);
	
	for (ic=1; ic < nc; ic++) {
	    sol[ic]=0.0;
	}

	for (ic=0; ic < nc-1; ic++) {
	    for (id=0; id <= ic; id++) {
		sol[id+1] -= roots[ic]*sol[id];
	    }
	}

	if (NULL != coef) sf_floatwrite(sol,nc,coef);

	for (i1=0; i1 < n; i1++) {
	    dat[i1] = 0.;
	    x = o+i1*d;
	    xp = 1.0;
	    for (ic=0; ic < nc; ic++) {
		dat[i1] += xp*sol[ic];
		xp *= x;
	    }
	}

	sf_floatwrite(dat,n,out);
    }

    exit(0);
}
