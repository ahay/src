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
#include "freqlet.h"

int main(int argc, char* argv[])
{
    int i1, n1, i2, n2, niter, iter, rect=1;
    bool adj, verb, *known;
    const char *type, *oper;
    sf_complex *z0, *x, *y, *y2=NULL;
    float *m=NULL, perc, eps;
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
	case 't':
	    coper = freqlet_lop;

	    if (NULL == (type = sf_getstring("type"))) type="haar";
	    /* [haar,linear,biortho] type of the seislet transform */

	    freqlet_init(n1, true, true, type[0]);
	    break;
	case 's':
	    coper = sin_smooth;

	    if (!sf_getint("rect",&rect)) rect=1;
	    /* smoothing radius (for oper=s) */
	    break;
	default:
	    sf_error("Unknown operator \"%s\"",oper);
    }


    x = sf_complexalloc(n1);
    y = sf_complexalloc(n1);
    known = sf_boolalloc(n1);

    if (niter > 0) { /* do inversion */
	if (NULL != sf_getstring ("mask")) {
	    /* missing data interpolation */
	    
	    mask = sf_input("mask");
	    if (SF_FLOAT != sf_gettype(mask)) sf_error("Need float mask");
	    m = sf_floatalloc(n1);
	 
	    switch (oper[0]) {
		case 't':
		    y2 = sf_complexalloc(n1);

		    if (!sf_getfloat("perc",&perc)) perc=50.;
		    /* percentage for thresholding (used when oper=t and niter > 0) */

		    sf_sharpen_init(n1,perc);
		    break;
		case 'd':
		    for (i1=0; i1 < n1; i1++) {
			x[i1] = sf_cmplx(0.,0.);
		    }
		    break;
		case 'c':
		    sf_mask_init(known);
		    break;
		case 's':
		    sf_mask_init(known);
		    y2 = sf_complexalloc(n1);

		    if (!sf_getfloat("eps",&eps)) eps=1./n1;
		    /* scaling for shaping inversion */

		    sf_cconjgrad_init(n1,n1,n1,n1, 
				      eps,
				      1.e-7,
				      true,false);
		    break;
	    }
	}
    }

    for (i2=0; i2 < n2; i2++) {
	switch(oper[0]) {
	    case 'd':
	    case 'c':
		sin_init(z0[i2]);
		break;
	    case 't':
		freqlet_setz(z0[i2]);
		break;
	    case 's':
		sinsmooth_init(z0[i2],n1,rect);
		break;
	}

	if (niter > 0) {
	    if (NULL != mask) { /* missing data interpolation */
		sf_floatread(m,n1,mask);
		for (i1=0; i1 < n1; i1++) {
		    known[i1] = (bool) (m[i1] != 0.);
		}
		sf_complexread(y,n1,in);

		switch (oper[0]) {
		    case 'd':
			sf_csolver(coper, sf_ccgstep, n1, n1, y, x, niter,
				   "known", known, "x0", y, "verb", verb, "end");
			break;
		    case 'c':	
			sf_csolver_prec(sf_cmask_lop, sf_ccgstep, coper, 
					n1, n1, n1, y, y, niter, 0.,
					"verb", verb,"end");
			break;
		    case 's':
			sf_cconjgrad(NULL, sf_cmask_lop, coper, 
				     y2, y, y, niter);
			break;
		    case 't':
			for (iter=0; iter < niter; iter++) {
			    coper(false,false,n1,n1,x,y2);
			    for (i1=0; i1 < n1; i1++) {
				if (known[i1]) y2[i1]=y[i1];
			    }
			    coper(true,false,n1,n1,x,y2);
			    sf_csharpen(x);
			    sf_cweight_apply(n1,x);
			}
			coper(false,false,n1,n1,x,y);
			break;
		}
	    } else {
		sf_error("Unknown operation");
	    }
	} else {
	    sf_complexread(x,n1,in);

	    if (adj) {
		coper(true,false,n1,n1,y,x);
	    } else {
		coper(false,false,n1,n1,x,y);
	    }
	}

	sf_complexwrite(y,n1,out);
    }


    exit(0);
}
    
