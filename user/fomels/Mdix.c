/* Convert RMS to interval velocity using LS and shaping regularization. */
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

#include "smoothder.h"

int main(int argc, char* argv[])
{
    int i, niter, nd, dim, n1, n2, i1, i2;
    int n[SF_MAX_DIM], box[SF_MAX_DIM];
    float **vr, **vi, **wt, **v0, wti;
    char key[6];
    sf_file vrms, vint, weight, vout;

    sf_init(argc,argv);
    vrms = sf_input("in");
    vint = sf_output("out");
    
    if (NULL != sf_getstring("weight")) {
	weight = sf_input("weight");
    } else {
	weight = NULL;
    }

    dim = sf_filedims (vrms,n);

    nd = 1;
    for (i=0; i < dim; i++) {
	nd *= n[i];
    }
    n1 = n[0];
    n2 = nd/n1;

    for (i=0; i < dim; i++) { 	 
	snprintf(key,6,"rect%d",i+1); 	 	 
	if (!sf_getint(key,box+i)) box[i]=1; 	 
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
    } 	 
	 
    smoothder_init(nd, dim, box, n);

    vr = sf_floatalloc2(n1,n2);
    vi = sf_floatalloc2(n1,n2);
    wt = sf_floatalloc2(n1,n2);
    v0 = sf_floatalloc2(n1,n2);

    sf_floatread(vr[0],nd,vrms);

    if (NULL != weight) {
	sf_floatread(wt[0],nd,weight);
    } else {
	for (i1=0; i1 < nd; i1++) {
	    wt[0][i1] = 1.0f;
	}
    }

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
    }
    
    sf_repeat_lop(false,true,nd,nd,v0[0],vr[0]);
    smoothder(niter, wt[0], vr[0], vi[0]);
 
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

/* 	$Id$	 */
