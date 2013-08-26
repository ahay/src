/* Gardner's DMO for regularly sampled 2-D data (slow method) 

   The input/ouput is (time,offset,midpoint).
*/
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
    int it, nt, ih, nh, ix, nx, ib, nb, id, nd, *fold, CDPtype;
    float dt, dh, dx, t0, h0, x0, t, h, x, eps, db, b0, b, sinb, cosb;
    float ***slice, ***tstr, ***hstr, ***xstr, ***slice2, *sum, sample;
    sf_file in, out;
    
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nx)) sf_error("No n3= in input");

    nd = nt*nh*nx;

    CDPtype=1;
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dx)) sf_error("No d3= in input");
    CDPtype=0.5+dh/dx;
    if (0 == CDPtype) CDPtype=1;
    if (1 != CDPtype) sf_warning("CDPtype=%d",CDPtype);
    sf_putint(out,"CDPtype",1);

    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&x0)) sf_error("No o3= in input");

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    slice  = sf_floatalloc3(nt,nh,nx);
    tstr   = sf_floatalloc3(nt,nh,nx);
    hstr   = sf_floatalloc3(nt,nh,nx);
    xstr   = sf_floatalloc3(nt,nh,nx); 
    slice2 = sf_floatalloc3(nt,nh,nx);

    sum = sf_floatalloc(nd);
    fold = sf_intalloc(nd);

    if (!sf_getint("nb",&nb)) nb=171;   /* number of angles */
    if (!sf_getfloat("b0",&b0)) b0=-85; /* first angle */
    if (!sf_getfloat("db",&db)) db=1;   /* angle increment */

    b0 *= SF_PI/180;
    db *= SF_PI/180;

    warp3_init(nt, t0, dt,
	       nh, h0, dh,
               nx, x0, dx,
	       nt, nh, nx, eps); 

    for (id=0; id < nd; id++) {
	sum[id] = 0.0f;
	fold[id] = 0;
    }

    sf_floatread(slice[0][0],nd,in);

    for (ib=0; ib < nb; ib++) {
	sf_warning("angle %d of %d;",ib+1,nb);

	b = b0+ib*db;
	sinb = sinf(b);
	cosb = cosf(b);
	
	for (ix=0; ix < nx; ix++) {
	    x = x0+ix*dx;   
	    for (ih=0; ih < nh; ih++) {
		h = h0 + ih*dh + (dh/CDPtype)*(ix%CDPtype);
		for (it=0; it < nt; it++) {
		    t = t0 + it*dt;
		    tstr[ix][ih][it] = t*cosb;
		    hstr[ix][ih][it] = h*cosb;
		    xstr[ix][ih][it] = x-h*sinb;
		}
	    }
	}

	warp3(slice,tstr,hstr,xstr,slice2);

	for (id=0; id < nd; id++) {
	    sample = slice2[0][0][id];

	    if (sample != 0.0f) {
		sum[id] += sample;
		fold[id]++;
	    }
	}
    }
    sf_warning(".");

    for (id=0; id < nd; id++) {
	if (fold[id] > 0) sum[id] /= fold[id];
    }

    sf_floatwrite(sum,nd,out);
    exit(0);
}
