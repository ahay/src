/* Muting.

Takes: < cmp.rsf > muted.rsf

Data is smoothly weighted inside the mute zone.
The weight is zero for t <       (x-x0) * slope0
The weight is one  for t >  tp + (x-x0) * slopep

The signs are reversed for inner=y.
*/

#include <rsf.h>

#include "mutter.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i2,i3, CDPtype;
    bool abs, half, inner;
    float tp, slope0, slopep, o1,d1,o2,d2,d3, x,x0,x1, v0, *data;
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
    if (!sf_getfloat("v0",&v0)) v0=1.45; 
    if (!sf_getfloat("slope0",&slope0)) slope0=1./v0;
    if (!sf_getfloat("slopep",&slopep)) slopep=slope0;
    if (!sf_getfloat("x0",&x1)) x1=0.;

    if (!sf_getbool("abs",&abs)) abs=true;
    /* if y, use absolute value |x-x0| */

    if (!sf_getbool("inner",&inner)) inner=false;
    /* if y, do inner muter */

    data = sf_floatalloc(n1);

    mutter_init(n1,o1,d1,abs,inner);

    for (i3=0; i3 < n3; i3++) { 
	x0= o2 + (d2/CDPtype)*(i3%CDPtype) - x1;
	for (i2=0; i2 < n2; i2++) { 
	    x = x0+i2*d2;

	    sf_floatread (data,n1,in);
	    mutter (tp,slope0,slopep, x, data);
	    sf_floatwrite (data,n1,out);
	}
    }

    sf_close();
    exit(0);
}

/* 	$Id: Mmutter.c,v 1.7 2004/04/19 21:51:16 fomels Exp $	 */
