/* Invert a warping function.

Takes: < warp.rsf > invwarp.rsf
*/

#include <rsf.h>

#include "stretch.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, n2, i2, m2;
    float *inp, *outp;
    float o1, d1, o2, d2, eps;
    map str;
    sf_file in, out, other;

    sf_init (argc, argv);
    in = sf_input("in");
    other = sf_input("other");
    out = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    m2 = sf_leftsize(in,1);

    if(!sf_histint(other,"n1",&n2)) sf_error ("No n1= in other");
    if(!sf_histfloat(other,"d1",&d2)) sf_error ("No d1= in other");
    if(!sf_histfloat(other,"o1",&o2)) o2 = 0.;

    sf_putint(out,"n1",n2);
    sf_putfloat(out,"d1",d2);
    sf_putfloat(out,"o1",o2);

    if (!sf_getfloat("eps",&eps)) eps=0.01; 
    /* Smoothness parameter */

    inp = sf_floatalloc (n1);
    outp = sf_floatalloc (n2);

    str = stretch_init (n2, o2, d2, n1, eps, false);

    for (i2=0; i2 < m2; i2++) {
	sf_floatread (inp,n1,in);

	stretch_define (str, inp);

	for (i1=0; i1 < n1; i1++) {
	    inp[i1] = o1 + i1*d1;
	}
	
	stretch_apply (str, inp, outp);

/*
  for (i1=0; i1 < n2; i1++) {
  outp[i1] -= (o2 + i1*d2);
  }
*/

	sf_floatwrite(outp,n2,out);
    }


    exit (0);
}
    
/* 	$Id$	 */

