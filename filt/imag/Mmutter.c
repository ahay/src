
/*
  Data is weighted by sine squared inside a mute zone.
  The weight is zero above the line	t <       x * slope0
  The weight is one after the line     t >  tp + x * slopep
  Suggested defaults: slopep = slope0= 1./1.45 sec/km;  tp=.150 sec
  
  Defaults:
  slope0 = 1./1.45 = .69	sec/km
  slopep = slope0
  tp     = .150		sec
*/

#include <rsf.h>

#include "mutter.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i2,i3, CDPtype;
    float tp,tm, slope0, slopep, slopem, o1,d1,o2,d2, x, x0, *data;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);

    if (!sf_histint(in,"CDPtype",&CDPtype)) CDPtype=1;

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
    
    if (!sf_getfloat("tp",&tp)) tp=0.150;
    if (!sf_getfloat("tm",&tm)) tm=0.150;
    if (!sf_getfloat("slope0",&slope0)) slope0=1./1.45;
    if (!sf_getfloat("slopep",&slopep)) slopep=slope0;
    if (!sf_getfloat("slopem",&slopem)) slopem=0.; 

    data = sf_floatalloc(n1);

    mutter_init(n1,o1,d1);

    for (i3=0; i3 < n3; i3++) { 
	x0= o2 + (d2/CDPtype)*(i3%CDPtype);
	for (i2=0; i2 < n2; i2++) { 
	    x = x0+i2*d2;

	    sf_read (data,sizeof(float),n1,in);
	    mutter (tp,slope0,slopep, x, data);
	    sf_write (data,sizeof(float),n1,out);
	}
    }

    exit(0);
}

