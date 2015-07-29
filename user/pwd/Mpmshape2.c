/* Missing data interpolation in 2-D using plane-wave destruction and shaping regularization. */
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
#include "pwdsl.h"

int main(int argc, char* argv[])
{
    int i, niter, n1, n2, n12, i3, n3, rect1, rect2, order;
    float *mm, *dd, **pp, lam;
    bool *known;
    sf_file in, out, dip, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint("niter",&niter)) niter=100;
    /* number of iterations */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    pp = sf_floatalloc2(n1,n2);
    mm = sf_floatalloc(n12);
    known = sf_boolalloc(n12);
    
    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
	dd = sf_floatalloc(n12);
    } else {
	mask = NULL;
	dd = NULL;
    }

    if (!sf_getint("rect1",&rect1)) rect1=3;
    if (!sf_getint("rect2",&rect2)) rect2=3;
    /* smoothing radius */
    
    pwdsl_init(n1,n2,order,rect1,rect2,0.01);
    pwdsl_set(pp);
    sf_mask_init(known);
    
    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d",i3+1,n3);

	sf_floatread(mm,n12,in);

	if (NULL != mask) {
	    sf_floatread(dd,n12,mask);
	} else {
	    dd = mm;
	}

	/* figure out scaling and make known data mask */
	lam = 0.;
	for (i=0; i < n12; i++) {
	    if (dd[i] != 0.) {
		known[i] = true;
		lam += 1.;
	    } else {
		known[i] = false;
	    }
	}
	lam = sqrtf(lam/n12);

	/* read dip */
	sf_floatread(pp[0],n12,dip);

	sf_conjgrad_init(n12, n12, n12, n12, lam, 10*FLT_EPSILON, true, true); 
	sf_conjgrad(NULL,sf_mask_lop,pwdsl_lop,dd,mm,mm,niter);
	sf_conjgrad_close();

	sf_floatwrite (mm,n12,out);
    }

    exit(0);
}
