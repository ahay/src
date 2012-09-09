/* Inverse 3-D warping */
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
#include "warp3.h"

int main(int argc, char* argv[])
{
    bool inv;
    int nt, ny, nx, n1, n2, n3, i4, n4, ntxy;
    float o1, d1, o2, d2, o3, d3, eps;
    float ***slice, ***tstr, ***ystr, ***xstr, ***slice2;
    sf_file in, out, warp;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    warp = sf_input("warp");

   if (!sf_getbool("inv",&inv)) inv=true;
    /* inversion flag */

    if (inv) {
       if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
       if (!sf_histint(in,"n2",&ny)) sf_error("No n2= in input");
       if (!sf_histint(in,"n3",&nx)) sf_error("No n3= in input");
    } else {
       if (!sf_histint(warp,"n1",&nt)) sf_error("No n1= in input");
       if (!sf_histint(warp,"n2",&ny)) sf_error("No n2= in input");
       if (!sf_histint(warp,"n3",&nx)) sf_error("No n3= in input");

       if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
       if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
       if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    }
    n4 = sf_leftsize(in,3);

    ntxy= nt*ny*nx;

    if (inv && !sf_getint("n1",&n1)) n1=nt;
    if (!sf_getfloat("d1",&d1) && !sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_getfloat("o1",&o1) && !sf_histfloat(in,"o1",&o1)) o1=0.;
    
    if (inv && !sf_getint("n2",&n2)) n2=ny;
    if (!sf_getfloat("d2",&d2) && !sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_getfloat("o2",&o2) && !sf_histfloat(in,"o2",&o2)) o2=0.;

    if (inv && !sf_getint("n3",&n3)) n3=nx;
    /* output samples - for inv=y */
    if (!sf_getfloat("d3",&d3) && !sf_histfloat(in,"d3",&d3)) d3=1.;
    /*( d1=1 output sampling - for inv=y )*/
    if (!sf_getfloat("o3",&o3) && !sf_histfloat(in,"o3",&o3)) o3=0.;
    /*( o1=0 output origin - for inv=y )*/ 

    if (inv) {
       sf_putint(out,"n1",n1);
       sf_putfloat(out,"d1",d1);
       sf_putfloat(out,"o1",o1);

       sf_putint(out,"n2",n2);
       sf_putfloat(out,"d2",d2);
       sf_putfloat(out,"o2",o2);

       sf_putint(out,"n3",n3);
       sf_putfloat(out,"d3",d3);
       sf_putfloat(out,"o3",o3);
    } else {
       sf_putint(out,"n1",nt);
       sf_putint(out,"n2",ny);
       sf_putint(out,"n3",nx);
    }

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    slice  = sf_floatalloc3(nt,ny,nx);
    tstr   = sf_floatalloc3(nt,ny,nx);
    xstr   = sf_floatalloc3(nt,ny,nx);
    ystr   = sf_floatalloc3(nt,ny,nx); 
    slice2 = sf_floatalloc3(n1,n2,n3);

    warp3_init(n1, o1, d1,
	       n2, o2, d2,
               n3, o3, d3,
	       nt, ny, nx, eps); 

    for (i4=0; i4 < n4; i4++) {
	sf_floatread(tstr[0][0],ntxy,warp);
        sf_floatread(ystr[0][0],ntxy,warp);
	sf_floatread(xstr[0][0],ntxy,warp);

        if (inv) {
            sf_floatread(slice[0][0],ntxy,in);
	    warp3(slice,tstr,ystr,xstr,slice2);
	    sf_floatwrite (slice2[0][0],n1*n2*n3,out);
        } else {	    
            sf_floatread(slice2[0][0],n1*n2*n3,in);
	    fwarp3(slice2,tstr,ystr,xstr,slice); 
	    sf_floatwrite (slice[0][0],ntxy,out);
        }
    }

    exit(0);
}
