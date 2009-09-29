/* Compute intercept and gradient by least squares. */
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

int main(int argc, char* argv[])
{
    bool half;
    int it,ih,ix, nt,nh,nx, CDPtype;
    float dt, dh, t0, h0, h, dy, sh, sh2, det;
    float *trace=NULL, *a=NULL, *b=NULL, *hh=NULL, *sx=NULL, *sxh=NULL;
    sf_file cmp=NULL, avo=NULL, offset=NULL;

    sf_init (argc,argv);
    cmp = sf_input("in");
    avo = sf_output("out");

    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(cmp,"n2",&nh)) sf_error("No n2= in input");
    nx = sf_leftsize(cmp,2);

    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");

    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	h0 = dh = 0.;
	hh = sf_floatalloc(nh);
    } else {
	offset = NULL;
	if (!sf_histfloat(cmp,"o2",&h0)) sf_error("No o2= in input");
	if (!sf_histfloat(cmp,"d2",&dh)) sf_error("No d2= in input");
	
	if (!sf_getbool("half",&half)) half=true;
	/* if y, the second axis is half-offset instead of full offset */
	
	if (half) {
	    dh *= 2.;
	    h0 *= 2.;
	}

	hh = NULL;
    }

    sf_putfloat(avo,"n2",2);

    CDPtype=1;
    if (sf_histfloat(cmp,"d3",&dy)) {
	CDPtype=0.5+0.5*dh/dy;
	if (1 != CDPtype) sf_histint(cmp,"CDPtype",&CDPtype);
    }
    sf_warning("CDPtype=%d",CDPtype);

    trace = sf_floatalloc(nt);
    a = sf_floatalloc(nt);
    b = sf_floatalloc(nt);
    sx = sf_floatalloc(nt);
    sxh = sf_floatalloc(nt);

    for (ix=0; ix < nx; ix++) {
	if (NULL != offset) sf_floatread(hh,nh,offset);

	sh=0.;
	sh2=0.;
	for (it=0; it < nt; it++) {
	    sx[it] = 0.;
	    sxh[it] = 0.;
	}

	for (ih=0; ih < nh; ih++) {
	    h = (NULL != offset)? hh[ih]: 
		h0 + ih * dh + (dh/CDPtype)*(ix%CDPtype);
	    sh += h;
	    sh2 += h*h;

	    sf_floatread(trace,nt,cmp); 

	    for (it=0; it < nt; it++) {
		sx[it] += trace[it];
		sxh[it] += h*trace[it];
	    }
	}

	det = nh*sh2-sh*sh;

	for (it=0; it < nt; it++) {
	    a[it] = (sx[it]*sh2-sxh[it]*sh)/det;
	    b[it] = (sxh[it]*nh-sx[it]*sh)/det;
	}

	sf_floatwrite(a,nt,avo);
	sf_floatwrite(b,nt,avo);
    } /* x */

    exit(0);
}
