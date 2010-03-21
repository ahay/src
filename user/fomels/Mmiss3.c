/* Missing data interpolation (N-dimensional) using shaping regularization. */
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

int main(int argc, char* argv[])
{
    int niter, dim, n[SF_MAX_DIM], rect[SF_MAX_DIM], n12, i;
    float *mm, *kk, *pp, eps;
    bool *known, force;
    char key[6];
    sf_file in, out, mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getint("niter",&niter)) niter=100;
    /* Number of iterations */

    dim = sf_filedims (in,n);
    n12 = 1;
    for (i=0; i < dim; i++) {
	n12 *= n[i];
	if (n[i] > 1) {
	    snprintf(key,6,"rect%d",i+1);
	    if (!sf_getint(key,rect+i)) rect[i]=1;
	} else {
	    rect[i]=1;
	}
    }

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
 
    if (!sf_getbool("force",&force)) force=true;
    /* if y, keep known values */

    sf_mask_init(known);
    sf_trianglen_init(dim, rect, n);

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

    eps = 0.;
    for (i=0; i < n12; i++) {
	if (known[i]) eps += 1.;
    }
    eps = sqrtf(eps/n12);

    if (force) {
	for (i=0; i < n12; i++) {
	    if (known[i]) kk[i]=mm[i];
	}
    }

    sf_conjgrad_init(n12, n12, n12, n12, eps, FLT_EPSILON, true, false);

 
    sf_conjgrad(NULL, sf_mask_lop, sf_trianglen_lop, pp, mm, mm, niter);

    if (force) {
	for (i=0; i < n12; i++) {
	    if (known[i]) mm[i]=kk[i];
	}
    }
    
    sf_floatwrite(mm,n12,out);

    exit (0);
}

/* 	$Id$	 */
