#include <rsf.h>

#include "eno2.h"

int main(int argc, char* argv[]) 
{
    int n1, n2, i1, i2, order;
    float **dat, ***der;
    eno2 ent;
    sf_file in, out, deriv;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    deriv = sf_output("der");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

    sf_putint(deriv,"n1",2);
    sf_putint(deriv,"n2",n1);
    sf_putint(deriv,"n3",n2);

    dat = sf_floatalloc2(n1,n2);
    der = sf_floatalloc3(2,n1,n2);

    sf_read(dat[0],sizeof(float),n1*n2,in);

    if (!sf_getint("order",&order)) sf_error("Need order=");

    ent = eno2_init(order,n1,n2);
    eno2_set (ent, dat);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    eno2_apply (ent, i1, i2, 0., 0.,dat[i2]+i1,der[i2][i1],BOTH);
	}
    }

    sf_write(dat[0],sizeof(float),n1*n2,out);
    sf_write(der[0][0],sizeof(float),n1*n2*2,deriv);

    sf_close();
    exit(0);
}

/* 	$Id: Testeno2.c,v 1.2 2004/03/22 05:43:24 fomels Exp $	 */

