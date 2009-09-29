/* 2-D missing data interpolation. */
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
#include "gaussshape2.h"

int main(int argc, char* argv[])
{
    int niter, nliter, n1, n2, n3, i3, n12, i;
    float *mm=NULL, *kk=NULL, *pp=NULL, filt1, filt2, a[3], eps;
    bool *known=NULL, shape, force;
    sf_file in=NULL, out=NULL, mask=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */
    if (!sf_getint("nliter",&nliter)) nliter=1;
    /* Number of reweighting iterations */

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    mm = sf_floatalloc(n12);
    kk = sf_floatalloc(n12);
    pp = sf_floatalloc(n12);
    known = sf_boolalloc(n12);

    if (NULL != sf_getstring("mask")) {
	/* optional input mask file for known data */
	mask = sf_input("mask");
    } else {
	mask = NULL;
    }

    if (!sf_getfloat("filt1",&filt1)) filt1=3.;
    if (!sf_getfloat("filt2",&filt2)) filt2=filt1;
    /* smoothing radius */

    if (!sf_getfloat("a0",a))   a[0] = (filt1*filt1-1.)/12.;
    if (!sf_getfloat("b0",a+1)) a[1] = 0.;
    if (!sf_getfloat("c0",a+2)) a[2] = (filt2*filt2-1.)/12.;
    /* initial Gaussian shape parameters */

    if (!sf_getfloat("eps",&eps)) eps=0.0001;
    /* regularization parameter */

    if (!sf_getbool("shape",&shape)) shape=false;
    /* if y, estimate shaping */

    if (!sf_getbool("force",&force)) force=true;
    /* if y, keep known values */

    sf_mask_init(known);
    gaussshape2_init(n1,n2);
    sf_conjgrad_init(n12, n12, n12, n12, eps, 1.e-9, true, false);

    gaussshape2_set2(a);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(mm,n12,in);

	if (NULL != mask) {
	    sf_floatread(kk,n12,mask);
	    for (i=0; i < n12; i++) {
		known[i] = (bool) (kk[i] != 0.);
	    }
	} else {
	    for (i=0; i < n12; i++) {
		known[i] = (bool) (mm[i] != 0.);
	    }
	}

	if (force) {
	    for (i=0; i < n12; i++) {
		if (known[i]) kk[i]=mm[i];
	    }
	}
 
	sf_conjgrad(NULL, sf_mask_lop, sf_freqfilt2_lop, pp, mm, mm, niter);

	if (shape) {
	    gaussshape2_set(a, mm, 100, nliter);
	    sf_conjgrad(NULL, sf_mask_lop, sf_freqfilt2_lop, pp, mm, mm, niter);
	}

	if (force) {
	    for (i=0; i < n12; i++) {
		if (known[i]) mm[i]=kk[i];
	    }
	}

	sf_floatwrite(mm,n12,out);
    }

    exit (0);
}
