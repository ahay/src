/* Smoothing in 1-D by simple regularization.

Takes: < data.rsf > smooth.rsf
*/

#include <rsf.h>

#include "tridiagonal.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, n2, ir, nr;
    float *trace, eps,*diag, *offd;
    tris slv;
    sf_file in, smooth;

    sf_init (argc, argv);
    in = sf_input("in");
    smooth = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");

    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* smoothness parameter */
    eps = (eps*eps-1.)/12.;

    if (!sf_getint("repeat",&nr)) nr=1;
    /* repeat smoothing */

    n2 = sf_leftsize(in,1);

    trace = sf_floatalloc (n1);

    diag = sf_floatalloc (n1);
    offd = sf_floatalloc (n1);

    for (i1=0; i1 < n1; i1++) {
	diag[i1] = 1.+ 2.*eps;
	offd[i1] = -eps;
    }
    diag[0] = diag[n1-1] = 1.+eps;

    slv = tridiagonal_init(n1);
    tridiagonal_define (slv,diag,offd);

    for (i2=0; i2 < n2; i2++) {
	sf_read(trace,sizeof(float),n1,in);

	for (ir=0; ir < nr; ir++){
	    tridiagonal_solve(slv,trace);
	}

	sf_write(trace,sizeof(float),n1,smooth);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Msmoothreg.c,v 1.1 2004/04/08 01:57:36 fomels Exp $	 */
