/* Missing dara interpolation in 2-D by Laplacian regularization.

Takes: < inp.rsf [mask=mask.rsf] > filled.rsf
*/

#include <rsf.h>
#include "lapfill.h"

int main(int argc, char* argv[])
{
    int n1,n2, n12, n3, i1, i3, niter;
    float *map, *msk;
    bool *known, grad;
    sf_file in, out, mask;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    n12 = n1*n2;

    if (!sf_getint("niter",&niter)) niter=200;
    /* number of iterations */
    if (!sf_getbool("grad",&grad)) grad=false;
    /* if y, use gradient instead of laplacian */

    map = sf_floatalloc(n12);
    known = sf_boolalloc(n12);

    if (NULL != sf_getstring("mask")) {
	mask = sf_input("mask");
	msk =  sf_floatalloc(n12);
    } else {
	mask = NULL;
	msk = map;
    }

    lapfill_init (n1,n2,grad);

    for (i3=0; i3 < n3; i3++) {
	sf_read(map,sizeof(float),n12,in);
	
	if (NULL != mask) sf_read(msk,sizeof(float),n12,mask);

	for (i1=0; i1 < n12; i1++) {
	    known[i1] = (msk[i1] != 0.);
	}

	lapfill(niter,map,known);

	sf_write(map,sizeof(float),n12,out);
    }

    sf_close();
    exit(0);
}
