/* Muting.
   
Data is smoothly weighted inside the mute zone.
The weight is zero for t <       (x-x0) * slope0
The weight is one  for t >  tp + (x-x0) * slopep

The signs are reversed for inner=y.
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
    int n1, n2, n3, i2,i3, CDPtype;
    bool abs, half, inner, hyper;
    float t0, tp, slope0, slopep, o1,d1,o2,d2,d3, x,x0,x1, v0, *data;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

    if (!sf_getbool("half",&half)) half=true;
    /* if y, the second axis is half-offset instead of full offset */

    if (half) {
	d2 *= 2.;
	o2 *= 2.;
    }

    CDPtype=1;
    if (sf_histfloat(in,"d3",&d3)) {
	CDPtype=0.5+0.5*d2/d3;
	if (CDPtype < 1) CDPtype=1;
	if (1 != CDPtype) sf_histint(in,"CDPtype",&CDPtype);
    } 	    
    sf_warning("CDPtype=%d",CDPtype);
    
    if (!sf_getfloat("tp",&tp)) tp=0.150;
    if (!sf_getfloat("t0",&t0)) t0=0.;
    if (!sf_getfloat("v0",&v0)) v0=1.45; 
    if (!sf_getfloat("slope0",&slope0)) slope0=1./v0;
    if (!sf_getfloat("slopep",&slopep)) slopep=slope0;
    if (!sf_getfloat("x0",&x1)) x1=0.;

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
	x0= o2 + (d2/CDPtype)*(i3%CDPtype) - x1;
	for (i2=0; i2 < n2; i2++) { 
	    x = x0+i2*d2;
	    if (hyper) x *= x;
	    
	    sf_floatread (data,n1,in);
	    mutter (tp,slope0,slopep, x, data);
	    sf_floatwrite (data,n1,out);
	}
    }
    
    exit(0);
}

/* 	$Id$	 */
