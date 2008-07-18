/* Missing data interpolation in 2-D using seislet transform. */
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
#include "seislet.h"

int main(int argc, char* argv[])
{
    int i, niter, nw, n1, n2, n12, i3, n3; 
    float *mm, *dd, **pp, eps;
    char *type;
    bool verb, *known, prec;
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
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (NULL == (type=sf_getstring("type"))) type="biorthogonal";
    /* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal  */

    if (!sf_getbool("prec",&prec)) prec = false;
    /* if y, apply preconditioning */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */

    pp = sf_floatalloc2(n1,n2);

    mm = sf_floatalloc(n12);
    dd = sf_floatalloc(n12);
    known = sf_boolalloc(n12);

    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }

    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }

    if (prec) {
        /* initialize inverse seislet transform */

	sf_mask_init(known);
    } else {
	seislet_init(n1,n2,false,true,eps,type[0]);
	seislet_set(pp);
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
	
	if (prec) {
//	    sf_solver_prec(sf_mask_lop, sf_cgstep, predict_lop, 
//			   n12, n12, n12, 
//			   mm, mm, niter, 0.,"verb", verb,"end");
	} else {
	    sf_solver(seislet_lop, sf_cgstep, n12, n12, mm, dd, niter,
		      "known", known, "x0", mm, "verb", verb, "end");
	}
	
	sf_cgstep_close();
	
	sf_floatwrite (mm,n12,out);
    }
    
    
    exit(0);
}

/* 	$Id: Mseismis2.c 	 */
