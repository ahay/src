/* Smoothing in 2-D by simple regularization.

Takes: < data.rsf > smooth.rsf
*/
#include <math.h>

#include <rsf.h>

#include "copy.h"
#include "igrad2.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, n2, n12, ir, nr;
    float *trace, eps, *out;
    sf_file in, smooth;

    sf_init (argc, argv);
    in = sf_input("in");
    smooth = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");

    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* smoothness parameter */
    eps = sqrtf((eps*eps-1.)/12.);

    if (!sf_getint("repeat",&nr)) nr=1;
    /* repeat smoothing */

    n2 = sf_leftsize(in,1);
    n12 = n1*n2;

    trace = sf_floatalloc (n12);
    out = sf_floatalloc (n12);

    sf_floatread(trace,n12,in);

    igrad2_init(n1,n2);

    for (ir=0; ir < nr; ir++){
	sf_solver_reg (copy_lop, sf_cgstep, igrad2_lop, 
		       2*n12, n12, n12, out, trace, 
		       2*n12, eps, "end");
	sf_cgstep_close();

	for (i1=0; i1 < n12; i1++) {
	    trace[i1] = out[i1];
	}
    }

    sf_floatwrite(trace,n12,smooth);

    sf_close();
    exit (0);
}

/* 	$Id: Msmoothreg2.c,v 1.1 2004/05/18 11:41:27 fomels Exp $	 */
