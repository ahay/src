/* T squared stretch of the time axis (for use in velocity continuation).

Takes: < input.rsf > output.rsf
*/

#include <math.h>
#include <float.h>

#include <rsf.h>

#include "fint1.h"

int main(int argc, char* argv[])
{
    fint1 str;
    bool inv;
    int i1,i2, n1,n2, n, dens, it, nw;
    float d1, o1, d2, o2, *trace, *stretched, t;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse stretching */
    if (!sf_getint("dens",&dens)) dens=1;
    /* axis stretching factor */

    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (o1 < FLT_EPSILON) o1=FLT_EPSILON;

    if (!sf_histfloat(in,"d1",&d1))  sf_error("No d1= in input");

    if (!inv) {
	if (!sf_getint("nout",&n)) n=dens*n1;
	/* output axis length (if inv=n) */
	sf_putint(out,"nin",n1);

	o2 = o1*o1;
	d2 = o1+(n1-1)*d1;
	d2 = (d2*d2 - o2)/(n-1);
    } else {
	if (!sf_histint(in,"nin",&n)) n=n1/dens;

	o2 = sqrtf(o1);
	d2 = (sqrtf(o1+(n1-1)*d1) - o2)/(n-1);
    }

    sf_putint(out,"n1",n);
    sf_putfloat(out,"o1",o2);
    sf_putfloat(out,"d1",d2);

    trace = sf_floatalloc(n1);
    stretched = sf_floatalloc(n);

    if (!sf_getint("extend",&nw)) nw=4;
    /* trace extension */
    str = fint1_init(nw,n1);

    for (i2=0; i2 < n2; i2++) {
	sf_read (trace,sizeof(float),n1,in);
	fint1_set(str,trace);

	for (i1=0; i1 < n; i1++) {
	    t = o2+i1*d2;
	    t = inv? t*t: sqrtf(t);
	    t = (t-o1)/d1;
	    it = t;
	    if (it >= 0 && it < n1) {
		stretched[i1] = fint1_apply(str,it,t-it,false);
	    } else {
		stretched[i1] = 0.;
	    }
	}

        sf_write (stretched,sizeof(float),n,out);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mt2stretch.c,v 1.3 2004/03/22 05:43:25 fomels Exp $	 */

