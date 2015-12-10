/* Convert RMS to interval velocity using LS and shaping regularization.

Takes: rect1= rect2= ...

rectN defines the size of the smoothing stencil in N-th dimension.
*/
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

#include <math.h>

#include <rsf.h>

#include "smoothshape.h"

int main(int argc, char* argv[])
{
    int niter, nd, n1, n2, i1, i2, rect1, rect2, order;
    float **vr, **vi, **wt, **v0, **dp, wti, lam;
    sf_file vrms, vint, weight, vout, dip;

    sf_init(argc,argv);
    vrms = sf_input("in");
    vint = sf_output("out");
    weight = sf_input("weight");
    dip = sf_input("dip");

    if (!sf_histint(vrms,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(vrms,"n2",&n2)) sf_error("No n2= in input");
    nd = n1*n2;
    if (sf_leftsize(vrms,2) > 1) sf_error("Can handle 2D only");

    if (!sf_getint("rect1",&rect1)) rect1=3;
    if (!sf_getint("rect2",&rect2)) rect2=3;
    /* smoothing radius */

    if (!sf_getfloat("lam",&lam)) lam=1.;
    /* operator scaling for inversion */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    vr = sf_floatalloc2(n1,n2);
    vi = sf_floatalloc2(n1,n2);
    wt = sf_floatalloc2(n1,n2);
    v0 = sf_floatalloc2(n1,n2);
    dp = sf_floatalloc2(n1,n2);

    sf_floatread(vr[0],nd,vrms);
    sf_floatread(wt[0],nd,weight);
    sf_floatread(dp[0],nd,dip);

    smoothshape_init(n1, n2, order, rect1, rect2, lam, dp);

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */

    wti = 0.;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    wti += wt[i2][i1]*wt[i2][i1];
	}
    }
    if (wti > 0.) wti = sqrtf(n1*n2/wti);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vr[i2][i1] *= vr[i2][i1]*(i1+1.); /* vrms^2*t - data */
	    wt[i2][i1] *= wti/(i1+1.); /* decrease weight with time */	 
	    v0[i2][i1] = -vr[i2][0];
	}
    }
    
    sf_repeat_lop(false,true,nd,nd,v0[0],vr[0]);
    smoothshape(niter, wt[0], vr[0], vi[0]);
 
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vi[i2][i1] -= v0[i2][i1];
	}
    }

    sf_repeat_lop(false,false,nd,nd,vi[0],vr[0]);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vr[i2][i1] = sqrtf(fabsf(vr[i2][i1]/(i1+1.0f)));
	    vi[i2][i1] = sqrtf(fabsf(vi[i2][i1]));
	}
    }

    sf_floatwrite(vi[0],nd,vint);

    if (NULL != sf_getstring("vrmsout")) {
	/* optionally, output predicted vrms */
	vout = sf_output("vrmsout");

	sf_floatwrite(vr[0],nd,vout);
    }

    exit(0);
}

/* 	$Id: Mdix.c 1131 2005-04-20 18:19:10Z fomels $	 */
