/* Missing data interpolation in 2-D using plane-wave destruction. */
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
#include "twoplane2.h"
#include "allp2.h"
#include "predict.h"
#include "predict2.h"

int main(int argc, char* argv[])
{
    int i, niter, nw, n1, n2, n12, np, i3, n3, nj1, nj2;
    float *mm, *dd, **pp, **qq;
    bool *known, verb, prec;
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

    if (!sf_getint("order",&nw)) nw=1;
    /* accuracy order */
    if (!sf_getint("nj1",&nj1)) nj1=1;
    /* antialiasing for first dip */
    if (!sf_getint("nj2",&nj2)) nj2=1;
    /* antialiasing for second dip */

    if (!sf_getbool("prec",&prec)) prec = false;
    /* if y, apply preconditioning */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    np = sf_leftsize(dip,2);

    pp = sf_floatalloc2(n1,n2);

    if (np > n3) {
	qq = sf_floatalloc2(n1,n2);
    } else {
	qq = NULL;
    }

    mm = sf_floatalloc(n12);
    dd = sf_floatalloc(n12);
    known = sf_boolalloc(n12);
    
    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }

    if (NULL != qq) {
	if (prec) {
	    predict2_init(n1,n2,0.0001,nw,pp,qq);
	    sf_mask_init(known);
	} else {
	    twoplane2_init(nw, nj1,nj2, n1,n2, pp, qq);
	}
    } else {
	if (prec) {
	    predict_init(n1,n2,0.0001,nw,1,false);
	    predict_set(pp);
	    sf_mask_init(known);
	} else {
	    allpass22_init(allpass2_init(nw, nj1, n1,n2, pp));
	}
    }

    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d",i3+1,n3);

	sf_floatread(mm,n12,in);

	if (NULL != mask) {
	    sf_floatread(dd,n12,mask);

	    for (i=0; i < n12; i++) {
		known[i] = (bool) (dd[i] != 0.);
		dd[i] = 0.;
	    }
	} else {
	    for (i=0; i < n12; i++) {
		known[i] = (bool) (mm[i] != 0.);
		dd[i] = 0.;
	    }
	}

	sf_floatread(pp[0],n12,dip);

	if (NULL != qq) {
	    sf_floatread(qq[0],n12,dip);

	    if (prec) {
		sf_solver_prec(sf_mask_lop, sf_cgstep, predict2_lop, 
			       n12, n12, n12, 
			       mm, mm, niter, 0.,"verb", verb,"end");
	    } else {
		sf_solver(twoplane2_lop, sf_cgstep, n12, n12, mm, dd, niter,
			  "known", known, "x0", mm, "verb", verb, "end");
	    }
	} else {
	    if (prec) {
		sf_solver_prec(sf_mask_lop, sf_cgstep, predict_lop, 
			       n12, n12, n12, 
			       mm, mm, niter, 0.,"verb", verb,"end");
	    } else {
		sf_solver(allpass21_lop, sf_cgstep, n12, n12, mm, dd, niter,
			  "known", known, "x0", mm, "verb", verb, "end");
	    }
	}
	sf_cgstep_close();

	sf_floatwrite (mm,n12,out);
    }

    exit(0);
}
