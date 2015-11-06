/* Convert RMS to interval velocity using LS and plane-wave construction. */
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

#include "smoothpwd.h"

int main(int argc, char* argv[])
{
    int ncycle, niter, n1, n2, n12, i1, i2, rect, order;
    float **vr, **vi, **wt, **v0, **p=NULL, wti, eps;
    bool verb;
    sf_file vrms, vint, weight, vout, slope;

    sf_init(argc,argv);
    vrms = sf_input("in");
    vint = sf_output("out");
    weight = sf_input("weight");

    if (!sf_histint(vrms,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(vrms,"n2",&n2)) sf_error("No n2= in input");
    n12 = n1*n2;

    slope = sf_input("slope");
    p = sf_floatalloc2(n1,n2);
    sf_floatread(p[0],n12,slope);
    sf_fileclose(slope);
    
    vr = sf_floatalloc2(n1,n2);
    vi = sf_floatalloc2(n1,n2);
    wt = sf_floatalloc2(n1,n2);
    v0 = sf_floatalloc2(n1,n2);

    sf_floatread(vr[0],n12,vrms);
    sf_floatread(wt[0],n12,weight);

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */

    if (!sf_getint("ncycle",&ncycle)) ncycle=10;
    /* number of cycles for anisotropic diffusion */

    if (!sf_getint("rect1",&rect)) rect=1;
    /* vertical smoothing radius */

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getfloat("eps",&eps)) eps=0.;
    /* regularization parameter */

    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    smoothpwd_init(n1,n2,0.0001,order,rect,p);
    
    wti = 0.;
    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    wti += wt[i2][i1]*wt[i2][i1];
	}
    }
    if (wti > 0.) wti = sqrtf(n1*n2/wti);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vr[i2][i1] *= vr[i2][i1]*(i1+1.0f);
	    wt[i2][i1] *= wti/(i1+1.0f); /* decrease weight with time */	 
	    v0[i2][i1] = -vr[i2][0];
	}
    }
    
    sf_repeat_lop(false,true,n12,n12,v0[0],vr[0]);

    smoothpwd(niter, ncycle, wt[0], vr[0], vi[0], verb, eps);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vi[i2][i1] -= v0[i2][i1];
	}
    }

    sf_repeat_lop(false,false,n12,n12,vi[0],vr[0]);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    vr[i2][i1] = sqrtf(fabsf(vr[i2][i1]/(i1+1.0f)));
	    vi[i2][i1] = sqrtf(fabsf(vi[i2][i1]));
	}
    }

    sf_floatwrite(vi[0],n12,vint);

    if (NULL != sf_getstring("vrmsout")) {
	/* optionally, output predicted vrms */
	vout = sf_output("vrmsout");

	sf_floatwrite(vr[0],n12,vout);
    }


    exit(0);
}

/* 	$Id: Mdix.c 803 2004-09-18 12:32:15Z fomels $	 */
