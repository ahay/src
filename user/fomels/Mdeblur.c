/* Non-stationary debluring by inversion */
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

#include "ntriangle1.h"

int main(int argc, char* argv[])
{
    int i1, i2, n1, n2, n12, **nr, **ns, nbox, niter, iter, nliter;
    float *data, *modl, *wght, eps;
    bool verb;
    sf_file in, out, rect;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    rect = sf_input("rect");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (SF_INT != sf_gettype(rect)) sf_error("Need int rect");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);
    n12 = n1*n2;

    data = sf_floatalloc(n12);
    modl = sf_floatalloc(n12);
    wght = sf_floatalloc(n12);
    nr = sf_intalloc2(n1,n2);
    ns = sf_intalloc2(n1,n2);

    sf_floatread(data,n12,in);
    sf_intread(nr[0],n12,rect);

    nbox=1;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    if (nbox < nr[i2][i1]) nbox = nr[i2][i1];
	    ns[i2][i1]=0;
	}
    }

    for (i1=0; i1 < n12; i1++) {
	wght[i1] = 1.;
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("nliter",&nliter)) nliter=1;
    /* number of nonlinear iterations */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=0.;
    /* regularization parameter */

    ntriangle1_init(nbox,n1,n2,nr,ns);
    sf_weight_init(wght);
    sf_hilbert_init(n1, 10, 1.);

    for (iter=0; iter < nliter; iter++) {
	sf_solver_prec(ntriangle1_lop,sf_cgstep,sf_weight_lop,
		       n12,n12,n12,modl,data,niter,eps,
		       "verb",verb,"end");
	sf_cgstep_close();

	for (i2=0; i2 < n2; i2++) {
	    sf_hilbert(modl+i2*n1,wght+i2*n1);
	    for (i1=0; i1 < n1; i1++) {
		wght[i1+i2*n1] = hypotf(modl[i1+i2*n1],wght[i1+i2*n1]);
	    }
	}
    }

    sf_floatwrite(modl,n12,out);

    exit(0);
}
