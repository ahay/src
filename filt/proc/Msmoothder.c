/* Smooth first derivative on the first axis.

Takes: < data.rsf > derivative.rsf
*/

#include <rsf.h>

#include "banded.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, n2;
    float *trace, *dtrace, d1, eps, *diag, **offd;
    bands slv;
    sf_file in, der;

    sf_init (argc, argv);
    in = sf_input("in");
    der = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) d1=1.;

    if (!sf_getfloat("eps",&eps)) eps=0.2;
    /* smoothness parameter */

    n2 = sf_leftsize(in,1);

    trace = sf_floatalloc (n1);
    dtrace = sf_floatalloc (n1);
    diag = sf_floatalloc (n1);
    offd = sf_floatalloc2 (n1,2);

    for (i1=0; i1 < n1; i1++) {
	diag[i1] = 1.+ 6.*eps;
	offd[0][i1] = -4.*eps;
	offd[1][i1] = eps;
    }
    diag[0] = diag[n1-1] = 1.+eps;
    diag[1] = diag[n1-2] = 1.+5.*eps;
    offd[0][0] = offd[0][n1-2] = -2.*eps;

    slv = banded_init(n1,2);
    banded_define (slv,diag,offd);

    for (i2=0; i2 < n2; i2++) {
	sf_read(trace,sizeof(float),n1,in);

	/* smooth */
	banded_solve(slv,trace);

	/* differentiate */
	for (i1=0; i1 < n1-1; i1++) {
	    dtrace[i1] = trace[i1+1]-trace[i1];
	    dtrace[i1] /= d1;
	}
	dtrace[n1-1] = dtrace[n1-2];

	sf_write(dtrace,sizeof(float),n1,der);
    }

    exit (0);
}

/* 	$Id: Msmoothder.c,v 1.3 2003/10/01 22:45:55 fomels Exp $	 */
