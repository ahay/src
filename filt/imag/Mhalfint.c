/* Half-order integration or differentiation.

Takes: < in.rsf > out.rsf
*/

#include <rsf.h>

#include "halfint.h"

int main(int argc, char* argv[])
{
    bool adj, inv;
    int n1,n2,i2;
    float rho, *pp;
    sf_file in, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getbool("adj",&adj)) adj=false;
    /* If y, apply adjoint */
    if (!sf_getbool("inv",&inv)) inv=false;
    /* If y, do differentiation instead of integration */
    if (!sf_getfloat("rho",&rho)) rho = 1.-1./n1;
    /* Leaky integration constant */

    if (adj) {
	sf_warning("%s half-order differentiation",inv? "anticausal":"causal");
    } else {
	sf_warning("%s half-order integration",inv? "anticausal":"causal");
    }

    pp = sf_floatalloc(n1);

    halfint_init (adj, inv, n1, rho);

    for (i2=0; i2 < n2; i2++) {
	sf_read (pp,sizeof(float),n1,in);
	halfint (pp);
	sf_write (pp,sizeof(float),n1,out);
    }

    exit(0);
}

