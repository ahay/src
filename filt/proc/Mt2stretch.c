/* T squared stretch of the time axis (for use in velocity continuation).

Takes: < input.rsf > output.rsf
*/

#include <math.h>
#include <float.h>

#include <rsf.h>

#include "stretch.h"

int main(int argc, char* argv[])
{
    map str;
    bool inv;
    int i1,i2, n1,n2, n, dens;
    float d1, o1, d2, o2, eps, *trace, *s, *stretched, t;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse stretching */
    if (!sf_getfloat("eps",&eps)) eps=0.01;
    /* smoothness parameter */
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
    s = sf_floatalloc(n1);

    for (i1=0; i1 < n1; i1++) {
	t = o1+i1*d1;
	s[i1] = inv? sqrtf(t): t*t;
    }

    str = stretch_init (n, o2, d2, n1, eps);
    stretch_define (str, s);

    for (i2=0; i2 < n2; i2++) {
	sf_read (trace,sizeof(float),n1,in);
        stretch_apply (str, trace, stretched);
        sf_write (stretched,sizeof(float),n,out);
    }

    exit (0);
}

/* 	$Id: Mt2stretch.c,v 1.1 2003/12/04 05:14:31 fomels Exp $	 */

