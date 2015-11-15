/* Muting.
   
Data is smoothly weighted inside the mute zone.
The weight is zero for t <       (x-x0) * slope0
The weight is one  for t >  tp + (x-x0) * slopep

The signs are reversed for inner=y.

July 2015 program of the month:
http://ahay.org/blog/2015/07/10/program-of-the-month-sfmutter/
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

#include "mutter.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i2,i3, CDPtype, nh2;
    bool abs, half, inner, hyper, nan;
    float t0, tp, slope0, slopep, o1,d1,o2,d2,d3, x,x1, v0, *data, *off;
    sf_file in, out, offset;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */
    
     if (!sf_getbool("nan",&nan)) nan=false;
    /* if y, put  nans instead of zeros */
    
    CDPtype=1;
    if (NULL != sf_getstring("offset")) {
	offset = sf_input("offset");
	nh2 = sf_filesize(offset);
	if (nh2 != n2 && nh2 != n2*n3) sf_error("Wrong dimensions in offset");

	off = sf_floatalloc(nh2);	
	sf_floatread (off,nh2,offset);
	sf_fileclose(offset);
    } else {
	if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
	if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

	if (sf_histfloat(in,"d3",&d3)) {
	    CDPtype=half? 0.5+d2/d3 : 0.5+0.5*d2/d3;
	    if (CDPtype < 1) {
		CDPtype=1;
	    } else if (1 != CDPtype) {
		sf_histint(in,"CDPtype",&CDPtype);
	    	sf_warning("CDPtype=%d",CDPtype);
	    }
	} 

	nh2 = n2;
	off = sf_floatalloc(nh2);
	for (i2 = 0; i2 < n2; i2++) {
	    off[i2] = o2 + i2*d2; 
	}

	offset = NULL;
    }
    
    if (!sf_getfloat("tp",&tp)) tp=0.150; /* end time */
    if (!sf_getfloat("t0",&t0)) t0=0.;    /* starting time */
    if (!sf_getfloat("v0",&v0)) v0=1.45;  /* velocity */
    if (!sf_getfloat("slope0",&slope0)) slope0=1./v0; /* slope */
    if (!sf_getfloat("slopep",&slopep)) slopep=slope0; /* end slope */
    if (!sf_getfloat("x0",&x1)) x1=0.; /* starting space */

    if (!sf_getbool("abs",&abs)) abs=true;
    /* if y, use absolute value |x-x0| */

    if (!sf_getbool("inner",&inner)) inner=false;
    /* if y, do inner muter */

    if (!sf_getbool("hyper",&hyper)) hyper=false;
    /* if y, do hyperbolic mute */

    if (hyper) {
	slope0 *= slope0;
	slopep *= slopep;
    }

    data = sf_floatalloc(n1);

    mutter_init(n1,o1-t0,d1,abs,inner,hyper);
    
    for (i3=0; i3 < n3; i3++) { 
	for (i2=0; i2 < n2; i2++) { 
	    x = (nh2 == n2)? off[i2] + (d2/CDPtype)*(i3%CDPtype) : 
		off[i3*n2+i2];	    
	    x -= x1;
	    if (half) x *= 2.;
	    if (hyper) x *= x;
	    
	    sf_floatread (data,n1,in);
	    mutter (tp,slope0,slopep, x, data, nan);
	    sf_floatwrite (data,n1,out);
	}
    }


    exit(0);
}

/* 	$Id$	 */
