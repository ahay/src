/* Taking first derivative along the fast axis using ENO interpolation.

Takes: < input.rsf > deriv.rsf

*/

#include <rsf.h>

#include "eno.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, n2, order;
    float *trace, f, f1, d1;
    eno ent;
    sf_file in, der;

    sf_init (argc, argv);
    in = sf_input("in");
    der = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) d1=1.;

    n2 = sf_leftsize(in,1);

    if(!sf_getint("accuracy",&order)) order=4;
    /* Interpolation accuracy order */

    trace = sf_floatalloc (n1);
    ent = eno_init (order, n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(trace,n1,in);

	eno_set (ent, trace); 
	for (i1=0; i1 < n1; i1++) {
	    eno_apply (ent, i1, 0., &f, &f1, DER);
	    trace[i1] = f1/d1; 
	}

	sf_floatwrite(trace,n1,der);
    }


    exit (0);
}

/* 	$Id$	 */

