/* Add a perturbation to the warping function.

Takes: < warp1.rsf add=dwarp.rsf > warp2.rsf
*/

#include <math.h>

#include <rsf.h> 

#include "eno.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, n2, m2, order, i;
    float *first, *second, o1, d1, x, f, f1, t;
    eno ent;
    sf_file in, sum, add;

    sf_init (argc, argv);
    in = sf_input("in");
    add = sf_input("add");
    sum = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) sf_error ("No d1= in input");
    if(!sf_histfloat(in,"o1",&o1)) o1 = 0.;

    if(!sf_histint(in,"n2",&m2)) m2 = 1;

    if(!sf_getint("accuracy",&order)) {
	/* Interpolation accuracy order */
	order = 2;
    } else if (order < 1 || order > 4) {
	sf_error ("accuracy must be between 1 and 4");
    }
    order *= 2;

    if (!sf_getint("m1",&n2)) n2=n1*2;
    /* Trace pading */

    ent = eno_init(order,n2);

    first = sf_floatalloc (n1); 
    second = sf_floatalloc (n2);

    for (i2=0; i2 < m2; i2++) {
	sf_read(first,sizeof(float),n1,add);
	sf_read(second,sizeof(float),n1,in);

	for (i1=n1; i1 < n2; i1++) {
	    second[i1] = (7.*second[i1-1]-5.*second[i1-2]+second[i1-3])/3.;
	}

	eno_set (ent,second);

	for (i1=0; i1 < n1; i1++) {
	    t = (o1+i1*d1);
	    x = (first[i1] + t - o1)/d1;
	    i = floorf(x); x -= i;
	    eno_apply (ent, i, x, &f, &f1, FUNC);
	    first[i1] += f;
	}

	sf_write(first,sizeof(float),n1,sum);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mwarpadd.c,v 1.4 2004/03/22 05:43:24 fomels Exp $	 */
