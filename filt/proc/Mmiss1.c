/* Missing data interpolation in 1-D.

Takes: < inp.rsf [filtin=filt.rsf] > out.rsf
*/

#include <rsf.h>

#include "mis1.h"

int main(int argc, char* argv[])
{
    int i1, i2, n1, n2, na;
    float *xx, *aa;
    bool *known;
    sf_file in, out, filt;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (NULL != sf_getstring("filtin")) {
	filt = sf_input("filtin");
	if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in filtin");

	aa = sf_floatalloc(na);
	sf_floatread(aa,na,filt);
    }  else {
	aa = sf_floatalloc(2);
	aa[0] = 1.;
	aa[1] = -1.;
    }

    xx = sf_floatalloc(n1);
    known = sf_boolalloc(n1);

    mis1_init(n1,na,aa);

    for (i2=0; i2 < n2; i2++) {
        sf_floatread (xx,n1,in);

	for (i1=0; i1 < n1; i1++) {
	    known[i1] = (xx[i1] != 0.);
	}

        mis1 (n1, xx, known);

	sf_floatwrite (xx,n1,out);
    }

    sf_close();
    exit(0);
}

