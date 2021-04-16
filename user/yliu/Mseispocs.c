/* Missing data interpolation using POCS added 2-D seislet transform. */
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
#include <rsfpwd.h>

int main(int argc, char* argv[])
{
    int i, niter, nw, n1, n2, n12, i1, i3, n3, iter, order; 
    float *mm, *dd, *dd2=NULL, **pp, *m=NULL, eps, perc,fact, *sknown;
    char *type;
    bool verb, *dknown;
    sf_file in, out, dip, mask=NULL, convex=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    dip = sf_input("dip");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;
    n3 = sf_leftsize(in,2);

    if (!sf_getint("niter",&niter)) niter=20;
    /* number of iterations */

    if (!sf_getint("order",&nw)) nw=1;
    /* [1,2,3] accuracy order */
    if (nw < 1 || nw > 3) 
	sf_error ("Unsupported nw=%d, choose between 1 and 3",nw);

    if (NULL == (type=sf_getstring("type"))) type="biorthogonal";
    /* [haar,linear,biorthogonal] wavelet type, the default is biorthogonal  */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* regularization parameter */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    pp = sf_floatalloc2(n1,n2);
    mm = sf_floatalloc(n12);
    dd = sf_floatalloc(n12);
    dknown = sf_boolalloc(n12);
    sknown = sf_floatalloc(n12);
    m = sf_floatalloc(n12);

    if (NULL != sf_getstring ("mask")) {
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }
    if (NULL != sf_getstring ("convex")) {
	convex = sf_input("convex");
    } else {
	convex = NULL;
	if (!sf_getfloat("perc",&perc)) perc=90.;
        /* percentage for smooth */
	if (!sf_getfloat("fact",&fact)) fact=0.5;
        /* percentage for smooth */
	
	sf_sharpen_init(n12,perc,fact);
    }

    if (NULL != mask) {
	sf_floatread(m,n12,mask);
	
	for (i=0; i < n12; i++) {
	    dknown[i] = (bool) (m[i] != 0.);
	}
    } else {
	for (i=0; i < n12; i++) {
	    dknown[i] = (bool) (dd[i] != 0.);
	}
    }
    if (NULL != convex) {
	sf_floatread(sknown,n12,convex);
    }  

    seislet_init(n1,n2,true,true,eps,order,type[0]);
    dd2 = sf_floatalloc(n12);

    seislet_set(pp);

    for (i3=0; i3 < n3; i3++) {
	sf_warning("slice %d of %d",i3+1,n3);

	sf_floatread(dd,n12,in);
	sf_floatread(pp[0],n12,dip);
	for (i1=0; i1 < n12; i1++) {
	    dd2[i1] = dd[i1];
	}

	/* seislet transform POCS */
	for (iter=0; iter < niter-1; iter++) {
	    if (verb)
		sf_warning("Iteration %d of %d",iter+1,niter);

	    seislet_lop(true,false,n12,n12,mm,dd2);
	    if (NULL != convex) {
		for (i1=0; i1 < n12; i1++) {
		    mm[i1]*=sknown[i1];
		}
	    } else {
		sf_sharpen(mm);
		sf_weight_apply(n12,mm);
	    }
	    seislet_lop(false,false,n12,n12,mm,dd2);
	    for (i1=0; i1 < n12; i1++) {
		if (dknown[i1]) dd2[i1]=dd[i1];
	    }
	}
	
	sf_floatwrite (dd2,n12,out);
    }

    exit(0);
}

/* 	$Id$	 */
