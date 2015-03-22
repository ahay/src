/* Missing data interpolation in 3-D using plane-wave shaping regularization. */
/*
  Copyright (C) 2015 University of Texas at Austin
  
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
#include "pwsmooth3.h"

int main(int argc, char* argv[])
{
    int i, niter, n1, n2, n3, n123, i4, n4, ns1,ns2, order1, order2;
    float *mm, *dd, *xx, ****pp, lam, eps;
    bool *known;
    sf_file in, out, dip, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    n123 = n1*n2*n3;
    n4 = sf_leftsize(in,3);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order1",&order1)) order1=1;
    if (!sf_getint("order2",&order2)) order2=1;
    /* accuracy order */

    pp = sf_floatalloc4(n1,n2,n3,2);
    mm = sf_floatalloc(n123);
    xx = sf_floatalloc(n123);
    known = sf_boolalloc(n123);
    
    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
	dd = sf_floatalloc(n123);
    } else {
	mask = NULL;
	dd = NULL;
    }

    if (!sf_getint("ns1",&ns1)) ns1=1;
    if (!sf_getint("ns2",&ns2)) ns2=1;
    /* smoothing radius */
    
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization */

    sf_mask_init(known);
    
    for (i4=0; i4 < n4; i4++) {
	sf_floatread(mm,n123,in);
	for (i=0; i < n123; i++) {
	    xx[i] = mm[i];
	}

	if (NULL != mask) {
	    sf_floatread(dd,n123,mask);
	} else {
	    dd = mm;
	}

	/* figure out scaling and make known data mask */
	lam = 0.;
	for (i=0; i < n123; i++) {
	    if (dd[i] != 0.) {
		known[i] = true;
		lam += 1.;
	    } else {
		known[i] = false;
	    }
	}
	lam = sqrtf(lam/n123);

	/* read dip */
	sf_floatread(pp[0][0][0],n123*2,dip);

	pwsmooth3_init(ns1,ns2, n1, n2, n3, order1, order2, eps, pp);

	sf_conjgrad_init(n123, n123, n123, n123, lam, 10*FLT_EPSILON, true, true); 
	sf_conjgrad(NULL,sf_mask_lop,pwsmooth3_lop,xx,mm,mm,niter);
	sf_conjgrad_close();

	pwsmooth3_close();

	sf_floatwrite (mm,n123,out);
    }

    exit(0);
}
