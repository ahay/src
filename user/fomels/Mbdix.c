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

#include "blockder.h"

int main(int argc, char* argv[])
{
    int niter, nd, n1, n2, i1, i2;
    float **vr, **vi, **wt, **v0, **bk, *tmp, wti, perc;
    sf_file vrms, vint, weight, vout, block;

    sf_init(argc,argv);
    vrms = sf_input("in");
    vint = sf_output("out");
    weight = sf_input("weight");
    block = sf_input("block");

    if (!sf_histint(vrms,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(vrms,1);
    nd = n1*n2;

    vr = sf_floatalloc2(n1,n2);
    vi = sf_floatalloc2(n1,n2);
    wt = sf_floatalloc2(n1,n2);
    v0 = sf_floatalloc2(n1,n2);
    bk = sf_floatalloc2(n1,n2);
    tmp = sf_floatalloc(n1);

    sf_floatread(vr[0],nd,vrms);
    sf_floatread(wt[0],nd,weight);
    sf_floatread(bk[0],nd,block);

    if (!sf_getfloat("perc",&perc)) perc=50.0;
    /* percentage for sharpening */

    blockder_init(n1, n2, perc, bk[0], wt[0]);

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
	    vr[i2][i1] *= vr[i2][i1]*(i1+1.0f); /* vrms^2*t - data */
	    wt[i2][i1] *= wti/(i1+1.0f); /* decrease weight with time */	 
	    v0[i2][i1] = -vr[i2][0];
	}
	sf_causint_lop(false,true,n1,n1,v0[i2],vr[i2]);
    }
    
    blockder(niter, vr[0], vi[0]);
 
    for (i2=0; i2 < n2; i2++) {
	sf_causint_lop(false,false,n1,n1,vi[i2],tmp);
	for (i1=0; i1 < n1; i1++) {
	    vi[i2][i1] = tmp[i1] - v0[i2][i1];
	}
	sf_causint_lop(false,false,n1,n1,vi[i2],vr[i2]);
    }

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
