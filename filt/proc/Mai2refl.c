/* 1-D convolution modeling.

Takes: < ai.rsf > reflectivity.rsf
*/
#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int nt, it, n2, i2;
    float imp1, imp2, *imp, *sig;
    sf_file ai, mod;

    sf_init (argc,argv);
    ai = sf_input("in");
    mod = sf_output("out");

    if (!sf_histint(ai,"n1",&nt)) sf_error("No n1= in input");
    n2 = sf_leftsize(ai,1);

    imp = sf_floatalloc (nt);
    sig = sf_floatalloc (nt);

    for (i2=0; i2 < n2; i2++) {
	for (it=0; it < nt; it++) {
	    sig[it]=0.;
	}

	sf_floatread(imp,nt,ai);

	imp1=imp[0];
	for (it=0; it < nt-1; it++) {
	    imp2 = imp[it+1];
	    sig[it] = (imp2-imp1)/(imp2+imp1);
	    imp1 = imp2;
	}
	sig[nt-1] = 0.;

	sf_floatwrite(sig,nt,mod);
    }

    exit (0);
}

/* 	$Id$	 */
