/* 1-D convolution modeling.

Takes: < ai.rsf > reflectivity.rsf
*/

#include <math.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int nt, it, n2, i2, is;
    float dt, f, imp1, imp2, a, rc;
    float *imp, *sig;
    sf_file ai, mod;

    sf_init (argc,argv);
    ai = sf_input("in");
    mod = sf_output("out");

    if (!sf_histint(ai,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(ai,"d1",&dt)) sf_error("No d1= in input");

    n2 = sf_leftsize(ai,1);

    if (!sf_getfloat ("f",&f)) f=40.;
    /* Peak frequency of Ricker wavelet */
    f *= SF_PI*dt;

    imp = sf_floatalloc (nt);
    sig = sf_floatalloc (nt);

    for (i2=0; i2 < n2; i2++) {
	for (it=0; it < nt; it++) {
	    sig[it]=0.;
	}

	sf_read(imp,sizeof(float),nt,ai);

	imp1=imp[0];
	for (it=0; it < nt-1; it++) {
	    imp2 = imp[it+1];
	    rc = (imp2-imp1)/(imp2+imp1);
	    for (is=0; is < nt; is++) {
		a = f*(is-it);
		a *= a;
		sig[is] += rc*(1.-2.*a)*expf(-a);
	    }
	    imp1 = imp2;
	}

	sf_write(sig,sizeof(float),nt,mod);
    }

    sf_close();
    exit (0);
}

/* 	$Id: Mai2refl.c,v 1.3 2004/03/22 05:43:24 fomels Exp $	 */
