/* 1-D spline interpolation

Takes: < input.rsf > interpolated.rsf
*/

#include <rsf.h>

#include "fint1.h"

int main(int argc, char* argv[])
{
    int nd, n1, i1, n2, i2, ix;
    float x, o1, d1, *table1=NULL, *trace, x0, dx;
    fint1 fnt;
    sf_file in, out, pattern;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&nd)) sf_error("Need n1= in input");

    if (!sf_histfloat(in,"d1",&dx)) sf_error("Need d1= in input");
    if (!sf_histfloat(in,"o1",&x0)) sf_error("Need o1= in input");
    n2 = sf_leftsize(in,1);
    
    if (NULL != sf_getstring("pattern")) {
	pattern = sf_input("pattern");
    } else {
	pattern = NULL;
    }

    if (!sf_getint("n1",&n1) && 
	(NULL== pattern ||
	 !sf_histint(pattern,"n1",&n1))) sf_error("Need n1=");
    /* Output grid size */
    if (!sf_getfloat("d1",&d1) && 
	(NULL== pattern ||
	 !sf_histfloat(pattern,"d1",&d1))) sf_error("Need d1=");
    /* Output sampling */
    if (!sf_getfloat("o1",&o1) &&
	(NULL== pattern ||
	 !sf_histfloat(pattern,"o1",&o1))) sf_error("Need o1=");
    /* Output origin */

    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);

    trace = sf_floatalloc(n1);
    table1 = sf_floatalloc(nd);

    fnt = fint1_init (4,nd);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(table1,nd,in);
	fint1_set(fnt,table1);
	    
	for (i1=0; i1 < n1; i1++) {
	    x = ((o1 + i1*d1)-x0)/dx;
	    ix = x;
	    trace[i1] = fint1_apply(fnt,ix,x-ix,false);
	}

	sf_floatwrite(trace,n1,out);
    }
    
    sf_close();
    exit(0);
}

/* 	$Id: Mfint1.c,v 1.1 2004/06/03 05:37:09 fomels Exp $	 */
