#include <rsf.h>

#include "causint.h"
#include "cgstep.h"
#include "bigsolver.h"

int main(int argc, char* argv[])
{ 
    int i1, n1, i2, n2, niter;
    float *trace, *dtrace, d1, eps;
    sf_file in, der;

    sf_init (argc, argv);
    in = sf_input("in");
    der = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error ("No n1= in input");
    if(!sf_histfloat(in,"d1",&d1)) d1=1.;

    if (!sf_getint("niter",&niter)) niter=100;
    if (!sf_getfloat("eps",&eps)) eps=0.2;
    eps *= n1;

    n2 = sf_leftsize(in,1);

    trace = sf_floatalloc (n1);
    dtrace = sf_floatalloc (n1);

    for (i2=0; i2 < n2; i2++) {
	sf_read(trace,sizeof(float),n1,in);

	solver_prec( causint_lop, cgstep, causint_lop, n1, n1, n1,
		     dtrace, trace, niter, eps, "verb", false, "end");
	cgstep_close();

	for (i1=0; i1 < n1; i1++) {
	    dtrace[i1] /= d1;
	}

	sf_write(dtrace,sizeof(float),n1,der);
    }

    exit (0);
}
