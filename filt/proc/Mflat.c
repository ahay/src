/* Moveout flattening. */
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

#include "dip2.h"
#include "predict.h"

int main (int argc, char *argv[])
{
    bool verb, sign, gauss;
    int n1,n2,n3, niter, nw, nj1, i3,i2,i1;
    float eps, lam, v0, d1, d2, o1, o2;
    float **u, **p, **v;
    sf_file in, out, dip;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;

    if (!sf_getint("niter",&niter)) niter=10;
    if (!sf_getfloat("eps",&eps)) eps=1.; 
    if (!sf_getfloat("lam",&lam)) lam=1.;

    eps = sqrtf(12*eps+1.);
    lam = sqrtf(12*lam+1.);
 
    if (!sf_getfloat("v0",&v0)) v0=0.;
    if (v0 != 0.) v0=d2/(v0*v0*d1);
    if (!sf_getint("nw",&nw)) nw=1;
    if (nw < 1 || nw > 3) sf_error ("nw must be between 1 and 3");
    if (!sf_getint("nj1",&nj1)) nj1=1;

    if(!sf_getbool("verb",&verb)) verb = false;
    if(!sf_getbool("sign",&sign)) sign = false;

    if (NULL != sf_getstring("dip")) {
	dip = sf_output("dip");
    } else {
	dip = NULL;
    }
    
    if (!sf_getbool("gauss",&gauss)) gauss = true;
    /* if y, use exact Gaussian for smoothing */

    /* initialize dip estimation */
    dip2_init(n1, n2, eps, lam, sign, gauss);
    predict_init (n1, n2, 0.0001, verb);

    u = sf_floatalloc2(n1,n2);
    p = sf_floatalloc2(n1,n2);
    v = sf_floatalloc2(n1,n2);

    for (i3=0; i3 < n3; i3++) {
	if (verb) fprintf(stderr,"cmp %d of %d\n",i3+1,n3);

	sf_floatread(u[0],n1*n2,in);

	for (i2=0; i2 < n2; i2++) {
	    /* initialize dip */
	    for(i1=0; i1 < n1; i1++) {
		p[i2][i1] = v0*(i2*d2+o2)/(i1*d1+o1+0.1*d1);
	    }
	}

	/* estimate dip */
	dip2(niter, nw, nj1, verb, u, p);
	predict_flat(u, v, p);

	sf_floatwrite(v[0],n1*n2,out);
	if (NULL != dip) sf_floatwrite(p[0],n1*n2,dip);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mflat.c,v 1.3 2004/06/25 18:08:42 fomels Exp $	 */
