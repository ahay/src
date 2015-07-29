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
#include "warp2.h"

int main(int argc, char* argv[])
{
    bool inv;
    int nt, nx, n1, n2, i3, n3, ntx;
    float o1, d1, o2, d2, eps, dt, dx, t0, x0;
    float **slice, **tstr, **xstr, **slice2;
    sf_file in, out, warp;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    warp = sf_input("warp");

   if (!sf_getbool("inv",&inv)) inv=true;
    /* inversion flag */

    if (inv) {
       if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
       if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    } else {
       if (!sf_histint(warp,"n1",&nt)) sf_error("No n1= in input");
       if (!sf_histint(warp,"n2",&nx)) sf_error("No n2= in input");

       if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
       if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    }
    n3 = sf_leftsize(in,2);

    ntx = nt*nx;

    if (inv) {
	if (!sf_getint("n1",&n1)) n1=nt;
	if (!sf_getfloat("d1",&d1) && !sf_histfloat(in,"d1",&d1)) d1=1.;
        /*( d1=1 output sampling - for inv=y )*/
	if (!sf_getfloat("o1",&o1) && !sf_histfloat(in,"o1",&o1)) o1=0.;
        /*( o1=0 output origin - for inv=y )*/ 

	if (!sf_getint("n2",&n2)) n2=nx;
	/* output samples - for inv=y */
	if (!sf_getfloat("d2",&d2) && !sf_histfloat(in,"d2",&d2)) d2=1.;
	/*( d2=1 output sampling - for inv=y )*/
	if (!sf_getfloat("o2",&o2) && !sf_histfloat(in,"o2",&o2)) o2=0.;
	/*( o2=0 output origin - for inv=y )*/ 

	sf_putint(out,"n1",n1);
	sf_putfloat(out,"d1",d1);
	sf_putfloat(out,"o1",o1);
	
	sf_putint(out,"n2",n2);
	sf_putfloat(out,"d2",d2);
	sf_putfloat(out,"o2",o2);
    } else {
	if (!sf_histfloat(in,"d1",&d1)) d1=1.;
	if (!sf_histfloat(in,"o1",&o1)) o1=0.;

	if (!sf_histfloat(in,"d2",&d2)) d2=1.;
	if (!sf_histfloat(in,"o2",&o2)) o2=0.;

	if (!sf_getfloat("d1",&dt) && !sf_histfloat(warp,"d1",&dt)) dt=d1;
	if (!sf_getfloat("o1",&t0) && !sf_histfloat(warp,"o1",&t0)) t0=o1;
	
	if (!sf_getfloat("d2",&dx) && !sf_histfloat(warp,"d2",&dx)) dx=d2;
	if (!sf_getfloat("o2",&x0) && !sf_histfloat(warp,"o2",&x0)) x0=o2;

	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"o1",t0);
	
	sf_putint(out,"n2",nx);
	sf_putfloat(out,"d2",dx);
	sf_putfloat(out,"o2",x0);
    }
    
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    slice  = sf_floatalloc2(nt,nx);
    tstr   = sf_floatalloc2(nt,nx);
    xstr   = sf_floatalloc2(nt,nx);

    slice2 = sf_floatalloc2(n1,n2);

    warp2_init(n1, o1, d1,
	       n2, o2, d2,
	       nt, nx, eps); 

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(tstr[0],ntx,warp);
	sf_floatread(xstr[0],ntx,warp);

        if (inv) {
            sf_floatread(slice[0],ntx,in);
	    warp2(slice,tstr,xstr,slice2);
	    sf_floatwrite (slice2[0],n1*n2,out);
        } else {
            sf_floatread(slice2[0],n1*n2,in);
	    fwarp2(slice2,tstr,xstr,slice);
	    sf_floatwrite (slice[0],ntx,out);
        }
    }

    exit(0);
}
