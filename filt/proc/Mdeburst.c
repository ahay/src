/* Remove bursty noise by IRLS. */

#include <rsf.h>

#include "irls.h"
#include "deburst.h"

int main(int argc, char* argv[])
{
    int n1, n2, i2, niter;
    char *norm;
    float *data, *model, eps;
    sf_file in, out;
    sf_weight weight=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    data = sf_floatalloc(n1);
    model = sf_floatalloc(n1);

    if (!sf_getint("niter",&niter)) niter=10;
    /* number of iterations */
    if (!sf_getfloat("eps",&eps)) eps=1.;
    /* regularization parameter */
    if (NULL == (norm = sf_getstring("norm"))) {
	/* norm to use in IRLS (cauchy,l1) */
	weight=cauchy;
    } else {
	sf_warning("got %s",norm);

	switch(norm[0]) {
	    case 'c': case 'C':
		weight=cauchy;
		break;
	    case 'l': case 'L':
		weight=l1;
		break;
	    default:
		sf_error("unknown norm %s",norm);
		break;
	}
    }

    irls_init(n1);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread (data,n1,in);

	deburst (n1, niter, weight, eps, data, model);

	sf_floatwrite (model,n1,out);
    }

    sf_close();
    exit(0);
}


