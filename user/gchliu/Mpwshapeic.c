/* Least Square Imaging condition using structure-based shaping regularization.

Takes: rect1= rect2= ...

rectN defines the size of the smoothing stencil in N-th dimension.
*/
/*
  Copyright (C) 2009 China University of Petroleum-Beijing
            and University of Texas at Austin
  
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
#include <rsfpwd.h>

#include "smoothdown.h"

int main(int argc, char* argv[])
{
    int niter, nd, n1, n2, rect1, rect2, i1, i2, order; 
    float **uw, **rf, **dw, **dp, lam, **wt;  
    sf_file upgw, refl, dwgw, dip, weight; 

    sf_init(argc,argv);
    upgw = sf_input("in");
    refl = sf_output("out");
    dwgw = sf_input("down");
    weight = sf_input("weight");
    dip = sf_input("dip");

    if (!sf_histint(upgw,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(upgw,"n2",&n2)) sf_error("No n2= in input");
    nd = n1*n2;

    if (!sf_getint("rect1",&rect1)) rect1=3;
    if (!sf_getint("rect2",&rect2)) rect2=3;
    /* smoothing radius */

    if (!sf_getfloat("lam",&lam)) lam=1.;
    /* operator scaling for inversion */
     
    if (!sf_getint("order",&order)) order=1;
    /* accuracy order */

    uw = sf_floatalloc2(n1,n2);
    dw = sf_floatalloc2(n1,n2);
    rf = sf_floatalloc2(n1,n2);
    wt = sf_floatalloc2(n1,n2);
    dp = sf_floatalloc2(n1,n2);
   

    sf_floatread(uw[0],nd,upgw);
    sf_floatread(dw[0],nd,dwgw); 
    
    sf_floatread(wt[0],nd,weight);
    sf_floatread(dp[0],nd,dip);
    /* sf_warning("lllllllllllll"); */
    for (i2=0; i2 < n2; i2++) {
    	for (i1=0; i1 < n1; i1++) { 
	    rf[i2][i1] = uw[i2][i1];
	}
    }

    smoothdown_init(n1, n2, order, rect1, rect2, lam, dp, dw[0]);

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */
	
    smoothdown(niter, wt[0], uw[0], rf[0]);
    sf_floatwrite(rf[0],nd,refl);

    exit(0);
}

/* 	$Id: Mdixshape.c 1131 2005-04-20 18:19:10Z fomels $	 */

