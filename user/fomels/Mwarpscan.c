/* Multicomponent data registration analysis. */
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

#include <string.h>
#include <math.h>
#include <float.h>

#include <rsf.h> 

#include "warpscan.h"

int main(int argc, char* argv[])
{ 
    bool shift, verb, cheb, sign;
    int n1, m[4], ntr, n2, order, ng, rect[4], niter, n2g, dim;
    float **inp, **oth, o1, d1, o2, d2, g0, dg, g1, o, d, *rat1;
    sf_file in, warped, other;

    sf_init (argc, argv);
    in = sf_input("in");
    warped = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    if(!sf_histint(in,"n2",&m[2])) m[2] = 1;
    if(!sf_histint(in,"n3",&m[3])) m[3] = 1;
    ntr = m[2]*m[3];

    if (m[3] > 1) {
	dim = 4;
    } else if (m[2] > 1) {
	dim = 3;
    } else {
	dim = 2;
    }

    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity flag */
    if (!sf_getbool("cheb",&cheb)) cheb=false;
    /* use Chebyshev scan */
    if (!sf_getbool("sign",&sign)) sign=false;
    /* use signed similarity */
    if (!sf_getint("ng",&ng)) ng=1;
    /* number of gamma values */
    if (!sf_getfloat("g0",&g0)) sf_error("Need g0=");
    /* gamma origin */
    if (!sf_getfloat("dg",&dg)) dg=g0;
    /* gamma sampling */
    g1=g0+(ng-1)*dg;
    if (cheb) {
	sf_putfloat(warped,"cmin",g0);
	sf_putfloat(warped,"cmax",g1);
    }

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    /* vertical smoothing */
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    /* gamma smoothing */
    if (!sf_getint("rect3",&rect[2])) rect[2]=1;
    /* in-line smoothing */
    if (!sf_getint("rect4",&rect[3])) rect[3]=1;
    /* cross-line smoothing */
    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */
    if (!sf_getbool("shift",&shift)) shift=false;
    /* use shift instead of stretch */

    other = sf_input("other");

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;
    sf_putint(warped,"n1",n2);
    sf_putfloat(warped,"d1",d2);
    sf_putfloat(warped,"o1",o2);

    m[0] = n2;
    m[1] = ng;

    sf_putint  (warped,"n2",ng);
    if (cheb) {
	sf_putfloat(warped,"d2",180.0/(ng-1));
	sf_putfloat(warped,"o2",0);
    } else {
	sf_putfloat(warped,"d2",dg);
	sf_putfloat(warped,"o2",g0);
    }

    if (dim > 3) {
	if(!sf_histfloat(in,"d3",&d)) d=1.;
	if(!sf_histfloat(in,"o3",&o)) o=0.;
	
	sf_putint  (warped,"n4",m[3]);
	sf_putfloat(warped,"d4",d);
	sf_putfloat(warped,"o4",o);
    } 

    if (dim > 2) {
	if(!sf_histfloat(in,"d2",&d)) d=1.;
	if(!sf_histfloat(in,"o2",&o)) o=0.;
	
	sf_putint  (warped,"n3",m[2]);
	sf_putfloat(warped,"d3",d);
	sf_putfloat(warped,"o3",o);
    }

    if(!sf_getint("accuracy",&order)) {
	/* [1-4] interpolation accuracy */
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;
    
    n2g = n2*ng*ntr;
    
    inp = sf_floatalloc2 (n1,ntr);
    oth = sf_floatalloc2 (n2,ntr);
    rat1 = sf_floatalloc (n2g);

    warpscan_init(n1,o1,d1,n2,o2,d2,ng,g0,dg,ntr,order,dim,m,rect,niter,shift,verb);
    
    sf_floatread(inp[0],n1*ntr,in);
    sf_floatread(oth[0],n2*ntr,other);

    warpscan(cheb,sign,inp,oth,rat1);

    sf_floatwrite(rat1,n2g,warped);

    exit (0);
}

/* 	$Id: Mwarpscan.c 744 2004-08-17 18:46:07Z fomels $	 */
