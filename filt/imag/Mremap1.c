/* 1-D ENO interpolation.

Takes: < in.rsf > out.rsf

*/

#include <rsf.h>

#include "eno.h"

int main(int argc, char* argv[])
{
    int n1, n2, i1, nn1, i2, order, i;
    float o1, d1, oo1, dd1, *tin, *tout, f, f1;
    eno map;
    sf_file in, out, pattern;

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");    
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getint("n1",&nn1)) nn1=n1;
    /* Number of output samples */
    if (!sf_getfloat("d1",&dd1)) dd1=d1;
    /* Output sampling */
    if (!sf_getfloat("o1",&oo1)) oo1=o1;
    /* Output origin */

    if (NULL != sf_getstring("pattern")) {
	pattern = sf_input("pattern");
    } else {
	pattern = NULL;
    }

    if (!sf_getint("n1",&nn1) && 
	(NULL== pattern ||
	 !sf_histint(pattern,"n1",&nn1))) nn1=n1;
    /* Output grid size */
    if (!sf_getfloat("d1",&dd1) && 
	(NULL== pattern ||
	 !sf_histfloat(pattern,"d1",&dd1))) dd1=d1;
    /* Output sampling */
    if (!sf_getfloat("o1",&oo1) &&
	(NULL== pattern ||
	 !sf_histfloat(pattern,"o1",&oo1))) oo1=o1;
    /* Output origin */

    sf_putint(out,"n1",nn1);
    sf_putfloat(out,"d1",dd1);
    sf_putfloat(out,"o1",oo1);

    tin = sf_floatalloc(n1);
    tout = sf_floatalloc(nn1);

    if (!sf_getint("order",&order)) order=3;
    /* Interpolation order */
    if (order > n1) order=n1;

    map = eno_init (order,n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(tin,n1,in);
	eno_set (map,tin);

	for (i1=0; i1 < nn1; i1++) {
	    f = (oo1+i1*dd1-o1)/d1; i=f; f -= i;
	    eno_apply(map, i, f, tout+i1, &f1, FUNC);
	}
	sf_floatwrite(tout,nn1,out);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mremap1.c,v 1.5 2004/04/19 21:51:16 fomels Exp $	 */
