/* Remove random shot gathers from a 2-D dataset.

Takes: < data.rsf > output.rsf mask=mask.rsf
*/
#include <float.h>

#include <rsf.h>

#include "randn.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i1, i2, i3, is, **known;
    float *trace, *zero, *chance, perc, rand, sum;
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    mask = sf_output("mask");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");

    sf_putint(mask,"n1",n2);
    sf_putint(mask,"n2",n3);
    sf_putint(mask,"n3",1);
    sf_setformat(mask,"native_int");

    if (!sf_getfloat("perc",&perc)) perc=0.75;
    /* how many shots to remove */

    trace = sf_floatalloc(n1);
    zero = sf_floatalloc(n1);
    known = sf_intalloc2(n2,n3);
    chance = sf_floatalloc(n2+n3);

    for (i1=0; i1 < n1; i1++) {
	zero[i1] = 0.;
    }

    random0 (n2+n3,chance);

    for (i3=0; i3 < n3; i3++) { /* half-offset */
	for (i2=0; i2 < n2; i2++) { /* midpoint */
	    is = i2 + i3; /* shot */
	    rand = chance[is];

	    sf_read (trace,sizeof(float),n1,in);
	    sum = 0.;
	    for (i1=0; i1 < n1; i1++) {
		sum += trace[i1]*trace[i1];
	    }

	    if (rand > perc && sum > FLT_EPSILON) {
		sf_write (trace,sizeof(float),n1,out);
		known[i3][i2] = 1;
	    } else {
		sf_write (zero,sizeof(float),n1,out);
		known[i3][i2] = 0;
	    }
	}
    }
    sf_write (known[0],sizeof(int),n2*n3,mask);

    exit(0);
}

