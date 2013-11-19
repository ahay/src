/* Wavenumber-domain Gardner's DMO for regularly sampled 2-D data 

   The input/ouput is (time,offset,midpoint wavenumber).
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
#include "warp2.h"

int main(int argc, char* argv[])
{
    bool inv;
    int it, nt, ih, nh, ik, nk, ib, nb, **fold;
    float dt, dh, dk, t0, h0, k0, t, h, k, eps, db, b, sinb, cosb, tanb;
    float **slice, **tstr, **hstr, **slice2, **sum, amp, sample;
    sf_file in, out;
    
    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_getbool("inv",&inv)) inv=true;
    /* inversion flag */

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nh)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&nk)) sf_error("No n3= in input");

    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&dh)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dk)) sf_error("No d3= in input");

    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&h0)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"o3",&k0)) sf_error("No o3= in input");

    dk *= 2*SF_PI;
    k0 *= 2*SF_PI;

    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* stretch regularization */

    slice  = sf_floatalloc2(nt,nh);
    tstr   = sf_floatalloc2(nt,nh);
    hstr   = sf_floatalloc2(nt,nh);
    slice2 = sf_floatalloc2(nt,nh);

    sum  = sf_floatalloc2(nt,nh);
    fold = sf_intalloc2(nt,nh);

    if (!sf_getint("nb",&nb)) nb=86;   /* number of angles */
    if (!sf_getfloat("db",&db)) db=1;   /* angle increment */
    db *= SF_PI/180;

    warp2_init(nt, t0, dt,
	       nh, h0, dh,
	       nt, nh, eps); 
	
    for (ik=0; ik < nk; ik++) {
	k = k0+ik*dk;   
	sf_warning("wavenumber %d of %d;",ik+1,nk);

	sf_floatread(slice[0],nt*nh,in);

	for (ih=0; ih < nh; ih++) {
	    for (it=0; it < nt; it++) {
		sum[ih][it] = slice[ih][it];
		fold[ih][it] = 1;
	    }
	}
	
	for (ib=1; ib < nb; ib++) {
	    b = ib*db;
	    sinb = sinf(b);
	    cosb = cosf(b);
	    tanb = sinb/cosb;

	    for (ih=0; ih < nh; ih++) {
		h = h0 + ih*dh;
		for (it=0; it < nt; it++) {
		    t = t0 + it*dt;
		    if (inv) {
			tstr[ih][it] = t*cosb;
			hstr[ih][it] = h*cosb;
		    } else {
			tstr[ih][it] = t/cosb;
			hstr[ih][it] = h/cosb;
		    }
		}
	    }

	    if (inv) {
		warp2(slice,tstr,hstr,slice2);
	    } else {	    
		fwarp2(slice,tstr,hstr,slice2); 
	    }

	    for (ih=0; ih < nh; ih++) {
		h = h0 + ih*dh;
		if (inv) {	
		    amp = 2.0f*cosf(k*h*sinb);
		} else {
		    amp = 2.0f*cosf(k*h*tanb);
		}

		for (it=0; it < nt; it++) {
		    sample = slice2[ih][it]*amp;

		    if (sample != 0.0f) {
			sum[ih][it]  += sample;
			fold[ih][it] += 2;
		    }
		}
	    }
	} /* b ended */
	
	for (ih=0; ih < nh; ih++) {
	    for (it=0; it < nt; it++) {
		if (fold[ih][it] > 1) sum[ih][it] /= fold[ih][it];
	    }
	}
    
	sf_floatwrite(sum[0],nt*nh,out);
    }
    sf_warning(".");

    exit(0);
}
