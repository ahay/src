/* Inverse 2-D warping */
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
#include "stretch4.h"

int main(int argc, char* argv[])
{
    map4 map1, map2;
    int nt, nx, n1, n2, i1, i2, i3, n3, ntx;
    float o1, d1, o2, d2, eps;
    float *trace1, *trace2;
    float **slice, **tstr, **xstr, **xstr1, **slice1, **slice2;
    sf_file in, out, warp;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    warp = sf_input("warp");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n1= in input");
    n3 = sf_leftsize(in,2);

    ntx = nt*nx;

    if (!sf_getint("n1",&n1)) n1=nt;
    if (!sf_getfloat("d1",&d1) && !sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_getfloat("o1",&o1) && !sf_histfloat(in,"o1",&o1)) o1=0.;

    if (!sf_getint("n2",&n2)) n2=nx;
    if (!sf_getfloat("d2",&d2) && !sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_getfloat("o2",&o2) && !sf_histfloat(in,"o2",&o2)) o2=0.;

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"o1",o1);

    sf_putint(out,"n2",n2);
    sf_putfloat(out,"d2",d2);
    sf_putfloat(out,"o2",o2);

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    slice  = sf_floatalloc2(nt,nx);
    tstr   = sf_floatalloc2(nt,nx);
    xstr   = sf_floatalloc2(nt,nx);

    trace1 = sf_floatalloc(n1);
    trace2 = sf_floatalloc(n2);

    xstr1  = sf_floatalloc2(nx,n1);
    slice1 = sf_floatalloc2(nx,n1);

    slice2 = sf_floatalloc2(n1,n2);

    map1 = stretch4_init (n1, o1, d1, nt, eps);
    map2 = stretch4_init (n2, o2, d2, nx, eps);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(slice[0],ntx,in);
	sf_floatread(tstr[0],ntx,warp);
	sf_floatread(xstr[0],ntx,warp);

	for (i2=0; i2 < nx; i2++) {
	    stretch4_define (map1,tstr[i2]);	    

	    stretch4_apply  (map1,slice[i2],trace1);	
	    for (i1=0; i1 < n1; i1++) {
		slice1[i1][i2] = trace1[i1];
	    }

	    stretch4_apply  (map1,xstr[i2],trace1);
	    for (i1=0; i1 < n1; i1++) {
		xstr1[i1][i2] = trace1[i1];
	    }
	}

	for (i1=0; i1 < n1; i1++) {
	    stretch4_define (map2,xstr1[i1]);
	    stretch4_apply  (map2,slice1[i1],trace2);
	    
	    for (i2=0; i2 < n2; i2++) {
		slice2[i2][i1] = trace2[i2];
	    }
	}

	sf_floatwrite (slice2[0],n1*n2,out);
    }

    exit(0);
}
