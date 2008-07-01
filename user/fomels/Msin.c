/* Simple operations with complex sinusoids */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

#include "sin.h"

int main(int argc, char* argv[])
{
    int i1, n1, i2, n2, niter;
    bool adj, verb, *known;
    char *oper;
    sf_complex *z0, *x, *y;
    float *m=NULL;
    sf_coperator coper=NULL;
    sf_file in, out, root, mask=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    root = sf_input("root");

    if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (sf_filesize(root) != n2) sf_error("Wrong dimensions in root");
    z0 = sf_complexalloc(n2);
    sf_complexread(z0,n2,root);
    sf_fileclose(root);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* adjoint flag */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getint("niter",&niter)) niter=0;
    /* number of iterations */

    if (NULL == (oper = sf_getstring("oper"))) oper="destruct";
    /* operation to perform */

    switch (oper[0]) {
	case 'd':
	    coper = sin_destruct;
	    break;
	case 'c':
	    coper = sin_construct;
	    break;
	default:
	    sf_error("Unknown operator \"%s\"",oper);
    }


    x = sf_complexalloc(n1);
    y = sf_complexalloc(n1);
    known = sf_boolalloc(n1);

    if (niter > 0) { /* do inversion */
	if (NULL != sf_getstring ("mask")) {
	    mask = sf_input("mask");
	    if (SF_FLOAT != sf_gettype(mask)) sf_error("Need float mask");
	    m = sf_floatalloc(n1);
	}

	for (i1=0; i1 < n1; i1++) {
	    y[i1] = 0.;
	}
    }

    for (i2=0; i2 < n2; i2++) {
	sf_complexread(x,n1,in);
	sin_init(z0[i2]);

	if (niter > 0) {
	    if (NULL != mask) {
		sf_floatread(m,n1,mask);
		for (i1=0; i1 < n1; i1++) {
		    known[i1] = (bool) (m[i1] != 0.);
		}
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    known[i1] = (bool) (crealf(x[i1]) != 0. && cimagf(x[i1]) != 0.);
		}
	    }

	    sf_csolver(coper, sf_ccgstep, n1, n1, x, y, niter,
		       "known", known, "x0", x, "verb", verb, "end");
	    sf_complexwrite(x,n1,out);
	} else {
	    if (adj) {
		coper(true,false,n1,n1,y,x);
	    } else {
		coper(false,false,n1,n1,x,y);
	    }

	    sf_complexwrite(y,n1,out);
	}
    }

    exit(0);
}
    
