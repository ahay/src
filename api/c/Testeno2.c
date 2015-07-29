#include "eno2.h"
#include "file.h"
#include "error.h"
#include "getpar.h"
#include "alloc.h"

int main(int argc, char* argv[]) 
{
    int n1, n2, i1, i2, order;
    float **dat, ***der;
    sf_eno2 ent;
    sf_file in=NULL, out=NULL, deriv=NULL;

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

    sf_floatread(dat[0],n1*n2,in);

    if (!sf_getint("order",&order)) sf_error("Need order=");

    ent = sf_eno2_init(order,n1,n2);
    sf_eno2_set (ent, dat);

    for (i2=0; i2 < n2; i2++) {
	for (i1=0; i1 < n1; i1++) {
	    sf_eno2_apply (ent, i1, i2, 0., 0.,dat[i2]+i1,der[i2][i1],BOTH);
	}
    }

    sf_floatwrite(dat[0],   n1*n2,out);
    sf_floatwrite(der[0][0],n1*n2*2,deriv);
    

    exit(0);
}
