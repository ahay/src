/* Automatic traveltime picking (yet another)

Takes: < input.rsf > output.rsf

*/

#include <rsf.h>

#include "pick0.h"

int main (int argc, char *argv[])
{
    int n1, n2, n3, i2, i3, order;
    float *trace;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    
    if (!sf_getint("order",&order)) order=4;
    /* Accuracy order */

    trace = sf_floatalloc (n1);

    pick0_init (n1, n2, order);

    for (i3=0; i3 < n3; i3++) {
	for (i2=0; i2 < n2; i2++) {
	    sf_floatread(trace,n1,in);
	    pick0_set (i2, trace);
	}

	for (i2=0; i2 < n2; i2++) {
	    pick0_delta (i2, trace);
	    sf_floatwrite(trace,n1,out);
	}
    }
    

    exit (0);
}

/* 	$Id$	 */
