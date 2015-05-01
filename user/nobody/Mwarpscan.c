/* Data registration analysis. */
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

#include "window1.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, m2, m3, n2, n, order, ng, ig, i0, w, nw, iw, jw, nw1, nw2;
    int wi;
    float *coord, **inp, *out, **oth, o1, d1, o2, d2, g0, dg, g;
    float *corr, *win1, *win2, a, b, a2, b2, ab, h, dw;
    bool taper, diff, verb, shift;
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

    if (!sf_getbool("shift",&shift)) shift=false;
    /* use shift instead of stretch */

    other = sf_input("other");

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;

    if (m3 > 1) {
	sf_putint  (warped,"n4",ng);
	sf_putfloat(warped,"d4",dg);
	sf_putfloat(warped,"o4",g0);
    } else if (m2 > 1) {
	sf_putint  (warped,"n3",ng);
	sf_putfloat(warped,"d3",dg);
	sf_putfloat(warped,"o3",g0);
    } else {
	sf_putint  (warped,"n2",ng);
	sf_putfloat(warped,"d2",dg);
	sf_putfloat(warped,"o2",g0);
    }

    m2 *= m3;
    n = n2*m2;

    if(!sf_getint("accuracy",&order)) {
	/* [1-4] interpolation accuracy */
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;

    if (!sf_getint("nw",&nw)) sf_error ("Need nw=");
    /* number of windows */
    if (!sf_getint("w",&w)) sf_error ("Need w=");
    /* window size */
    if (!sf_getfloat("h",&h)) h=0.5*(w-1);
    /* window overlap */
    if (!sf_getbool("taper",&taper)) taper=true;
    /* window tapering */
    if (!sf_getbool("verb",&verb)) verb=true;
    /* w */

    if (!sf_getbool("diff",&diff)) diff=false;
    /* if y, compute difference power insted of correlation */

    dw = (n2-w)/(nw-1.);
    nw1 = (0.5*w+1.)/dw;
    nw2 = (0.5*w-2.)/dw;
    nw += nw1 + nw2;

    sf_putint(warped,"n1",nw);
    sf_putfloat(warped,"o1",o2+(0.5*w+1.)*d2-nw2*dw*d2);
    sf_putfloat(warped,"d1",dw*d2);

    window1_init (h,dw);

    coord = sf_floatalloc (n2); 
    inp =   sf_floatalloc2 (n1,m2);
    out =   sf_floatalloc (n2);
    oth =   sf_floatalloc2 (n2,m2);

    corr =  sf_floatalloc (nw);

    win1 = sf_floatalloc (w);
    win2 = sf_floatalloc (w);
    win2 = sf_floatalloc (w);

    sf_prefilter_init (order, n1, order*10);     
    for (i2=0; i2 < m2; i2++) {
	sf_floatread(inp[i2],n1,in);
	sf_prefilter_apply (n1, inp[i2]);
    }
    sf_prefilter_close();

    sf_floatread(oth[0],m2*n2,other);
    sf_fileclose(other);

    for (ig=0; ig < ng; ig++) {
	if (verb) sf_warning("scanned %d of %d",ig+1,ng);

	g = g0 + ig*dg;

	for (i1=0; i1 < n2; i1++) {
	    coord[i1] = shift? o2+i1*d2+g: (o2+i1*d2)*g;
	}

	sf_int1_init (coord, o1, d1, n1, sf_spline_int, order, n2);

	for (i2=0; i2 < m2; i2++) {
	    sf_int1_lop (false,false,n1,n2,inp[i2],out);

	    for (iw=0; iw < nw; iw++) {
		a2=0.;
		b2=0.;
		ab=0.;

		if (iw < nw1) {
		    jw = 0;
		    wi = 0.5*(w+1.) + iw*dw;
		} else if (iw >= nw-nw2) {
		    jw = iw-nw1;
		    wi = 0.5*(w+1.) + (nw-1-iw)*dw; 
		} else {
		    jw = iw-nw1;
		    wi = w;
		}

		i0 = window1_apply(jw,wi,out,taper,taper,win1);
		i0 = window1_apply(jw,wi,oth[i2],taper,taper,win2);
		for (i1=0; i1 < wi; i1++) {
		    a = win1[i1];
		    b = win2[i1];
		    ab += a*b;
		    a2 += a*a;
		    b2 += b*b;
		}

		corr[iw] = diff? 
		    a2 + b2 - 2.*ab:
		    ab/sqrtf(a2*b2+FLT_EPSILON);
	    }

	    sf_floatwrite(corr,nw,warped);
	}
    }

    exit (0);
}

/* 	$Id$	 */
