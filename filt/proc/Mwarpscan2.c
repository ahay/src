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

#include "int1.h"
#include "interp_spline.h"
#include "prefilter.h"
#include "div2.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, m2, m3, n2, order, ng, ig, rect1, rect2, niter, n2g, i;
    float *coord, *inp, **out, o1, d1, o2, d2, g0, dg, g, dout, doth, o, d;
    float *rat1, *rat2, *num, *den, *oth;
    sf_file in, warped, other;

    sf_init (argc, argv);
    in = sf_input("in");
    warped = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    if(!sf_histint(in,"n2",&m2)) m2 = 1;
    if(!sf_histint(in,"n3",&m3)) m3 = 1;

    if (!sf_getint("ng",&ng)) ng=1;
    /* number of gamma values */
    if (!sf_getfloat("g0",&g0)) sf_error("Need g0=");
    /* gamma origin */
    if (!sf_getfloat("dg",&dg)) dg=g0;
    /* gamma sampling */
    if (!sf_getint("rect1",&rect1)) rect1=1;
    /* vertical smoothing */
    if (!sf_getint("rect2",&rect2)) rect2=1;
    /* gamma smoothing */
    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */

    other = sf_input("other");

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;
    sf_putint(warped,"n1",n2);
    sf_putfloat(warped,"d1",d2);
    sf_putfloat(warped,"o1",o2);

    sf_putint  (warped,"n2",ng);
    sf_putfloat(warped,"d2",dg);
    sf_putfloat(warped,"o2",g0);

    if (m3 > 1) {
	if(!sf_histfloat(in,"d3",&d)) d=1.;
	if(!sf_histfloat(in,"o3",&o)) o=0.;
	
	sf_putint  (warped,"n4",m3);
	sf_putfloat(warped,"d4",d);
	sf_putfloat(warped,"o4",o);
    } 

    if (m2 > 1 || m3 > 1) {
	if(!sf_histfloat(in,"d2",&d)) d=1.;
	if(!sf_histfloat(in,"o2",&o)) o=0.;
	
	sf_putint  (warped,"n3",m2);
	sf_putfloat(warped,"d3",d);
	sf_putfloat(warped,"o3",o);
    }

    m2 *= m3;

    if(!sf_getint("accuracy",&order)) {
	/* [1-4] interpolation accuracy */
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;

    coord = sf_floatalloc (n2); 
    inp =   sf_floatalloc (n1);
    out =   sf_floatalloc2 (n2,ng);
    oth =   sf_floatalloc (n2);

    n2g = n2*ng;

    rat1 = sf_floatalloc (n2g);
    rat2 = sf_floatalloc (n2g);
    num = sf_floatalloc (n2g);
    den = sf_floatalloc (n2g);

    prefilter_init (order, n1, order*10);     
    div2_init(n2, ng, (float) rect1,(float) rect2, niter, false);
    
    for (i2=0; i2 < m2; i2++) {
	sf_floatread(inp,n1,in);
	prefilter_apply (n1, inp);

	sf_floatread(oth,n2,other);

	doth = 0.;
	for (i1=0; i1 < n2; i1++) {
	    doth += oth[i1]*oth[i1];
	}
	doth = sqrtf(n2/doth);
	
	dout = 0.;
	for (ig=0; ig < ng; ig++) {
	    g = g0 + ig*dg;

	    for (i1=0; i1 < n2; i1++) {
		coord[i1] = (o2+i1*d2)*g;
	    }

	    int1_init (coord, o1, d1, n1, spline_int, order, n2);

	    int1_lop (false,false,n1,n2,inp,out[ig]);

	    for (i1=0; i1 < n2; i1++) {
		dout += out[ig][i1]*out[ig][i1];
	    }
	}
	dout = sqrtf(n2*ng/dout);

	for (ig=0; ig < ng; ig++) {
	    for (i1=0; i1 < n2; i1++) {
		i = ig*n2+i1;
		den[i] = out[ig][i1]*dout;
		num[i] = oth[i1]*dout;
	    }
	}
	div2(num,den,rat1);
	
	for (ig=0; ig < ng; ig++) {
	    for (i1=0; i1 < n2; i1++) {
		i = ig*n2+i1;
		num[i] = out[ig][i1]*doth;
		den[i] = oth[i1]*doth;
	    }
	}
	div2(num,den,rat2);
	
	for (i=0; i < n2g; i++) {
	    if (rat1[i] > 0.) {
		if (rat2[i] > 0. || -rat2[i] < rat1[i]) {
		    rat1[i] = sqrtf(fabsf(rat1[i]*rat2[i]));
		} else {
		    rat1[i] = -sqrtf(fabsf(rat1[i]*rat2[i]));
		}
	    } else {
		if (rat2[i] < 0. || rat2[i] < -rat1[i]) {
		    rat1[i] = -sqrtf(fabsf(rat1[i]*rat2[i]));
		} else {
		    rat1[i] = sqrtf(fabsf(rat1[i]*rat2[i]));
		}
	    }
	}
	
	sf_floatwrite(rat1,n2g,warped);
    }

    exit (0);
}

/* 	$Id: Mwarpscan.c 744 2004-08-17 18:46:07Z fomels $	 */
